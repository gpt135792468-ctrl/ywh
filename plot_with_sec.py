"""
完整 OTEC+SLTD 仿真 + SEC 分析
SEC (Specific Energy Consumption): 比能耗，单位产水所需的能量
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import fsolve
import warnings
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────
#  NH3 物性函数
# ─────────────────────────────────────────────
def NH3_psat(T):
    TK = T + 273.15
    return 10**(4.86886 - 1113.928/(TK - 10.409)) * 100

def NH3_rhol(T):
    return 638.6 - 1.353*T - 0.00310*T**2

def NH3_hf(T):
    return 200.0 + 4.476*T + 0.00520*T**2

def NH3_hg(T):
    return 1462.0 + 1.530*T - 0.00560*T**2

def NH3_sf(T):
    return 1.000 + 0.01583*T - 1.80e-5*T**2

def NH3_sg(T):
    return 5.616 - 0.00620*T - 5.0e-6*T**2

def NH3_isentropic_exp(s_in, T_out):
    sf = NH3_sf(T_out); sg = NH3_sg(T_out)
    hf = NH3_hf(T_out); hg = NH3_hg(T_out)
    x = np.clip((s_in - sf)/(sg - sf), 0, 1)
    return hf + x*(hg - hf), x

# ─────────────────────────────────────────────
#  海水物性函数
# ─────────────────────────────────────────────
def TEC(T):
    a = [0,-0.148759,-0.267408,1.080760,1.269056,-4.089591,-1.871251,7.438081,-3.536296]
    DT = sum(a[i]*(T/630)**(i) for i in range(9))
    T68 = T - DT
    val = 113586.363652 + 227272.7273*T68
    T48 = -113586.36365 + (max(val,0))**0.5
    return T68, T48

def cpw_func(S, T):
    T68, _ = TEC(T)
    t = T68 + 273.15
    sp = S * 1000
    A =  5.328 - 9.76e-2*sp + 4.04e-4*sp**2
    B = -6.913e-3 + 7.351e-4*sp - 3.15e-6*sp**2
    C =  9.6e-6 - 1.927e-6*sp + 8.23e-9*sp**2
    D =  2.5e-9 + 1.66e-9*sp - 7.125e-12*sp**2
    return 1000*(A + B*t + C*t**2 + D*t**3)

def BPEw(S, T):
    A = -4.584e-4*T**2 + 2.823e-1*T + 17.95
    B =  1.536e-4*T**2 + 5.267e-2*T +  6.56
    return A*S**2 + B*S

def hlat(S, T):
    hfg = 2.501e6 - 2.369e3*T + 2.678e-1*T**2 - 8.013e-3*T**3 - 2.079e-5*T**4
    return hfg*(1 - S/1000)

# ─────────────────────────────────────────────
#  SLTD 修复版（自适应松弛，正确结构）
# ─────────────────────────────────────────────
def sltd_fixed(N, Twwo, Tcwo, s, NEFe, NEFc, mww_in, mcw_in):
    n = N
    mf = np.zeros(n+1); md = np.zeros(n+1); S = np.zeros(n+1)
    Tf = np.zeros(n+1); Td = np.zeros(n+1)
    Te = np.zeros(n);   Tc = np.zeros(n)
    D  = np.zeros(n);   D2 = np.zeros(n)
    DT1= np.zeros(n);   DT2= np.zeros(n)
    BPE= np.zeros(n);   Tloss= np.zeros(n)
    cpe= np.zeros(n);   hfge = np.zeros(n)

    Tf[0]  = Twwo;  Td[n] = Tcwo
    mf[0]  = mww_in; md[n] = mcw_in
    S[0]   = s

    # 初始化 Te
    for i in range(n):
        Te[i] = Tf[0] - (i+1)*(Tf[0] - Td[n] + 2)/(n+1)

    ERR = 1.0; ite = 0
    while ERR > 1e-6*n and ite < 50000:
        omega = max(0.5 * 0.9995**ite, 0.05)  # 自适应松弛
        TE = Te.copy()

        # 正向：蒸发器
        for i in range(n):
            Tave = (Tf[i] + Te[i])/2
            BPE[i] = BPEw(S[i], Tave)
            dT_stage = max(Tf[i] - Te[i], 1e-6)
            DT1[i] = min(NEFe + BPE[i]/dT_stage, 0.99)
            cpe[i]  = cpw_func(S[i], Tave)
            hfge[i] = hlat(0, Tave)
            Tloss[i]= 0.15
            D[i] = max(mf[i]*cpe[i]*(1-DT1[i])*(Tf[i]-Te[i])/hfge[i], 0)
            mf[i+1] = max(mf[i] - D[i], 1e-6)
            S[i+1]  = S[i]*mf[i]/mf[i+1]
            Tf[i+1] = Te[i] + DT1[i]*(Tf[i] - Te[i])

        # 逆向：冷凝器
        for i in range(n):
            j = n-1-i
            Tc[j]  = Te[j] - Tloss[j]
            DT2[j] = NEFc
            Td[j]  = Tc[j] + DT2[j]*(Td[j+1] - Tc[j])
            Tave2  = (Td[j+1] + Tc[j])/2
            cpc_j  = cpw_func(0, Tave2)
            hfgc_j = hlat(0, Tave2)
            D2[j]  = md[j+1]*cpc_j*(1-DT2[j])*(Tc[j]-Td[j+1])/hfgc_j
            md[j]  = md[j+1] + D2[j]

        # 更新 Te（自适应松弛）
        ERR = 0; DT_total = 0
        for i in range(n):
            Te_new = Tf[i] - 0.5*(D[i]+D2[i])*hfge[i]/(cpe[i]*mf[i]*(1-DT1[i]))
            ERR += abs(1 - TE[i]/Te_new) if abs(Te_new) > 1e-12 else 0
            Te[i] = TE[i] + omega*(Te_new - TE[i])
            Tf[i+1] = Te[i] + DT1[i]*(Tf[i] - Te[i])
            DT_total += D[i]
        ite += 1

    return Tf[n], Td[0], DT_total, ite, ERR

# ─────────────────────────────────────────────
#  LMTD 残差
# ─────────────────────────────────────────────
def safe_lmtd(dT1, dT2):
    if dT1 <= 0 or dT2 <= 0:
        return -1e3
    if abs(dT1-dT2) < 1e-8:
        return dT1
    return (dT1-dT2)/np.log(dT1/dT2)

def residuals(x, Twwi, Tcwi, mww, mcw, mr, UAe, UAc, cpw_val, eta_t, eta_p):
    Teva,Tcon,Twwo,Tcwo = x
    h2 = NH3_hg(Teva); s2 = NH3_sg(Teva); h4 = NH3_hf(Tcon)
    pe = NH3_psat(Teva); pc = NH3_psat(Tcon)
    vl = 1/NH3_rhol(Tcon)
    wp = vl*(pe-pc)*1e3/eta_p
    h1 = h4 + wp/1000
    h3s,_ = NH3_isentropic_exp(s2, Tcon)
    h3 = h2 - eta_t*(h2-h3s)
    Qeva_r = mr*(h2-h1)*1e3; Qeva_w = mww*cpw_val*(Twwi-Twwo)
    Qcon_r = mr*(h3-h4)*1e3; Qcon_w = mcw*cpw_val*(Tcwo-Tcwi)
    LMTDe = safe_lmtd(Twwi-Teva, Twwo-Teva)
    LMTDc = safe_lmtd(Tcon-Tcwi, Tcon-Tcwo)
    return [Qeva_r-Qeva_w, Qeva_r-UAe*LMTDe, Qcon_r-Qcon_w, Qcon_r-UAc*LMTDc]

# ─────────────────────────────────────────────
#  主扫描
# ─────────────────────────────────────────────
Twwi=28; Tcwi=5; mcw=10; cpw_val=4020
UAe=50e3; UAc=50e3; eta_t=0.80; eta_p=0.75; eta_wp=0.80
mr=0.1; Hww=10; Hcw=25   # ★ Hcw=25m（液压扬程，非取水深度）
N_stage=4; S0=0.035; NEFe=0.02; NEFc=0.02
g=9.81

mww_vec = np.arange(10, 95, 5, dtype=float)
results = []

x0 = np.array([22., 10., 25., 7.])
for mww in mww_vec:
    sol = fsolve(residuals, x0, args=(Twwi,Tcwi,mww,mcw,mr,UAe,UAc,cpw_val,eta_t,eta_p),
                 full_output=True)
    xk, info, ier, msg = sol
    if ier != 1 or np.linalg.norm(info['fvec']) > 1e-3:
        results.append(None); continue

    Teva,Tcon,Twwo,Tcwo = xk
    if Teva>=Twwi or Tcon<=Tcwi or Twwo<=Teva or Tcwo>=Tcon:
        results.append(None); continue

    # Rankine 性能
    h2=NH3_hg(Teva); s2=NH3_sg(Teva); h4=NH3_hf(Tcon)
    pe=NH3_psat(Teva); pc=NH3_psat(Tcon)
    vl=1/NH3_rhol(Tcon); wp=vl*(pe-pc)*1e3/eta_p
    h1=h4+wp/1000
    h3s,_ = NH3_isentropic_exp(s2, Tcon)
    h3=h2-eta_t*(h2-h3s)
    Wt=mr*(h2-h3)*1e3; Wp=mr*wp
    Qeva=mr*(h2-h1)*1e3
    Wnet=Wt-Wp; eta_cycle=Wnet/Qeva

    Www_pump=mww*g*Hww/eta_wp
    Wcw_pump=mcw*g*Hcw/eta_wp
    Wnet_sys=Wnet-Www_pump-Wcw_pump

    # SLTD
    T1s,T2s,DT,ite,err = sltd_fixed(N_stage,Twwo,Tcwo,S0,NEFe,NEFc,mww,mcw)

    # ★ 正确 GOR：Q_total = SLTD 驱动热量 mww·cpw·(Twwo-Tcwo)
    Q_sltd = mww*cpw_val*(Twwo-Tcwo)
    hfg_local = 2501e3 - 2.369e3*Twwo
    GOR = DT*hfg_local/Q_sltd
    RR  = DT/mww   # 回收率

    # ─── SEC 比能耗 ───────────────────────────────────────────
    # SEC_thermal  = Q_sltd / DT        [J/kg] 热能比能耗
    # SEC_electric = W_pump_total / DT  [J/kg → 转 kWh/m³]
    #   W_pump_total = Www_pump + Wcw_pump（海水泵，Rankine发电抵消部分）
    #   若考虑 Wnet 作为内部供能：
    #     W_net_to_sltd = Wnet_sys（系统净功率，正值表示有余量）
    # 1 kWh/m³ = 3600e3 J / 1000 kg = 3600 J/kg (淡水密度≈1000 kg/m³)

    SEC_th  = Q_sltd / DT / 1e3            # [kJ/kg]
    SEC_el  = (Www_pump + Wcw_pump) / DT / 3600  # [kWh/m³]
    SEC_el_net = max(-Wnet_sys, 0) / DT / 3600   # 净电耗 [kWh/m³]（仅在Wnet_sys<0时有意义）

    results.append(dict(
        mww=mww,Teva=Teva,Tcon=Tcon,Twwo=Twwo,Tcwo=Tcwo,
        Wt=Wt,Wp=Wp,Wnet=Wnet,Wnet_sys=Wnet_sys,
        Qeva=Qeva,eta_cycle=eta_cycle,
        Www_pump=Www_pump,Wcw_pump=Wcw_pump,
        DT=DT,GOR=GOR,RR=RR,
        T1s=T1s,T2s=T2s,ite=ite,err=err,
        SEC_th=SEC_th,SEC_el=SEC_el,SEC_el_net=SEC_el_net,
        Q_sltd=Q_sltd
    ))
    x0 = xk.copy()

ok = [r for r in results if r is not None]
mww_ok  = np.array([r['mww']       for r in ok])
Teva_ok = np.array([r['Teva']      for r in ok])
Tcon_ok = np.array([r['Tcon']      for r in ok])
Twwo_ok = np.array([r['Twwo']      for r in ok])
Tcwo_ok = np.array([r['Tcwo']      for r in ok])
Wnet_ok = np.array([r['Wnet']      for r in ok])/1e3
Wsys_ok = np.array([r['Wnet_sys']  for r in ok])/1e3
Wt_ok   = np.array([r['Wt']        for r in ok])/1e3
eta_ok  = np.array([r['eta_cycle'] for r in ok])*100
DT_ok   = np.array([r['DT']        for r in ok])
GOR_ok  = np.array([r['GOR']       for r in ok])
RR_ok   = np.array([r['RR']        for r in ok])*100
Www_ok  = np.array([r['Www_pump']  for r in ok])/1e3
Wcw_ok  = np.array([r['Wcw_pump']  for r in ok])/1e3
Wp_ok   = np.array([r['Wp']        for r in ok])/1e3
T1s_ok  = np.array([r['T1s']       for r in ok])
Qsltd_ok= np.array([r['Q_sltd']    for r in ok])/1e3
SEC_th  = np.array([r['SEC_th']    for r in ok])
SEC_el  = np.array([r['SEC_el']    for r in ok])
SEC_el_net=np.array([r['SEC_el_net'] for r in ok])
ite_ok  = np.array([r['ite']       for r in ok])

# ─────────────────────────────────────────────
#  打印结果表
# ─────────────────────────────────────────────
print('='*90)
print(f"{'mww':>5} {'Twwo':>6} {'Tcwo':>6} {'Wnet':>7} {'Wsys':>7} "
      f"{'DT':>7} {'GOR':>7} {'RR%':>6} {'SEC_th':>9} {'SEC_el':>9} {'ite':>5}")
print(f"{'kg/s':>5} {'°C':>6} {'°C':>6} {'kW':>7} {'kW':>7} "
      f"{'kg/s':>7} {'-':>7} {'%':>6} {'kJ/kg':>9} {'kWh/m³':>9} {'-':>5}")
print('-'*90)
for r in ok:
    print(f"{r['mww']:5.0f} {r['Twwo']:6.2f} {r['Tcwo']:6.2f} "
          f"{r['Wnet']/1e3:7.3f} {r['Wnet_sys']/1e3:7.3f} "
          f"{r['DT']:7.4f} {r['GOR']:7.4f} {r['RR']*100:6.3f} "
          f"{r['SEC_th']:9.1f} {r['SEC_el']:9.4f} {r['ite']:5d}")
print('='*90)
print()
print("SEC 说明：")
print("  SEC_th   = Q_sltd / DT        [kJ/kg]   单位产水耗热能（越小越好）")
print("  SEC_el   = W_泵总 / DT         [kWh/m³]  单位产水耗电（海水泵）")
print("  Q_sltd   = mww·cpw·(Twwo-Tcwo) = SLTD 热驱动量")

# ─────────────────────────────────────────────
#  图 1：主性能 (2×3)
# ─────────────────────────────────────────────
C1=[0.00,0.45,0.74]; C2=[0.85,0.33,0.10]
C3=[0.47,0.67,0.19]; C4=[0.49,0.18,0.56]
C5=[0.93,0.69,0.13]
mk = dict(markersize=7, markerfacecolor='w', linewidth=2.0)

fig, axes = plt.subplots(2,3,figsize=(15,9))
fig.patch.set_facecolor('white')

ax=axes[0,0]
ax.plot(mww_ok,Wt_ok, '-o',color=C1,label='W_t (透平)',**mk)
ax.plot(mww_ok,Wnet_ok,'-s',color=C2,label='W_net (Rankine)',**mk)
ax.plot(mww_ok,Wsys_ok,'-d',color=C3,label='W_sys (含泵损)',**mk)
ax.axhline(0,color='k',lw=0.8,ls='--')
ax.set_xlabel('m_ww [kg/s]'); ax.set_ylabel('功率 [kW]')
ax.set_title('(a) 发电功率'); ax.legend(fontsize=8); ax.grid(True)

ax=axes[0,1]
ax.plot(mww_ok,DT_ok,'-o',color=C2,**mk)
ax.set_xlabel('m_ww [kg/s]'); ax.set_ylabel('D_T [kg/s]')
ax.set_title('(b) SLTD 淡水产量'); ax.grid(True)

ax=axes[0,2]
ax.plot(mww_ok,GOR_ok,'-^',color=C4,**mk)
ax.axhline(1,color='gray',lw=1,ls='--',label='GOR=1 参考线')
ax.set_xlabel('m_ww [kg/s]'); ax.set_ylabel('GOR [-]')
ax.set_title('(c) 造水比 GOR'); ax.legend(fontsize=9); ax.grid(True)

ax=axes[1,0]
ax.plot(mww_ok,RR_ok,'-v',color=C5,**mk)
ax.set_xlabel('m_ww [kg/s]'); ax.set_ylabel('回收率 RR [%]')
ax.set_title('(d) 产水回收率 (DT/mww)'); ax.grid(True)

ax=axes[1,1]
ax.plot(mww_ok,Teva_ok,'-o',color=C2,label='T_eva',**mk)
ax.plot(mww_ok,Tcon_ok,'-s',color=C1,label='T_con',**mk)
ax.plot(mww_ok,Twwo_ok,'--d',color=C2,label='T_wwo',lw=1.5,markersize=5)
ax.plot(mww_ok,Tcwo_ok,'--^',color=C1,label='T_cwo',lw=1.5,markersize=5)
ax.set_xlabel('m_ww [kg/s]'); ax.set_ylabel('温度 [°C]')
ax.set_title('(f) 关键温度'); ax.legend(fontsize=8,ncol=2); ax.grid(True)

ax=axes[1,2]
ax.plot(mww_ok,eta_ok,'-o',color=C1,**mk)
ax.set_xlabel('m_ww [kg/s]'); ax.set_ylabel('η [%]')
ax.set_title('(e) Rankine 循环热效率'); ax.grid(True)

fig.suptitle('闭式 OTEC + SLTD 系统 — 热侧流量敏感性', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('result_main.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ 图1保存：result_main.png")

# ─────────────────────────────────────────────
#  图 2：SEC 比能耗专项分析
# ─────────────────────────────────────────────
fig, axes = plt.subplots(2,2,figsize=(12,9))
fig.patch.set_facecolor('white')

# --- (a) SEC_th 热比能耗 ---
ax=axes[0,0]
ax.plot(mww_ok, SEC_th, '-o', color=C4, **mk)
# 参考线：理想单效蒸发 GOR=1 时 SEC_th = hfg ≈ 2430 kJ/kg
ax.axhline(2430, color='gray', ls='--', lw=1.2, label='GOR=1 理论极限 (≈2430 kJ/kg)')
ax.set_xlabel('m_ww [kg/s]'); ax.set_ylabel('SEC_th [kJ/kg]')
ax.set_title('(a) 热能比能耗 SEC_th\n= Q_sltd / D_T', fontsize=11)
ax.legend(fontsize=9); ax.grid(True)
ax.text(0.03, 0.95, 'SEC_th = hfg_local / GOR\n↓ mww增大 → GOR降 → SEC升',
        transform=ax.transAxes, fontsize=8, va='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# --- (b) SEC_el 电比能耗（海水泵） ---
ax=axes[0,1]
ax.plot(mww_ok, SEC_el, '-s', color=C2, label='SEC_el (泵总耗电)', **mk)
# 参考值：SWRO约3-5 kWh/m³，MED约1.5-2.5 kWh/m³
ax.axhline(1.5, color='green', ls='--', lw=1, label='MED 典型值 1.5 kWh/m³')
ax.axhline(3.0, color='orange', ls='--', lw=1, label='SWRO 典型值 3.0 kWh/m³')
ax.set_xlabel('m_ww [kg/s]'); ax.set_ylabel('SEC_el [kWh/m³]')
ax.set_title('(b) 电能比能耗 SEC_el\n= W_泵 / D_T', fontsize=11)
ax.legend(fontsize=8); ax.grid(True)

# --- (c) Q_sltd vs DT —— 热量–产水关系 ---
ax=axes[1,0]
ax2 = ax.twinx()
ax.plot(mww_ok, Qsltd_ok, '-o', color=C1, label='Q_sltd (SLTD驱动热)', **mk)
ax2.plot(mww_ok, DT_ok,   '-s', color=C2, label='D_T (产水量)', **mk)
ax.set_xlabel('m_ww [kg/s]')
ax.set_ylabel('Q_sltd [kW]', color=C1)
ax2.set_ylabel('D_T [kg/s]', color=C2)
ax.tick_params(axis='y', labelcolor=C1)
ax2.tick_params(axis='y', labelcolor=C2)
ax.set_title('(c) SLTD 热量输入 & 产水量', fontsize=11); ax.grid(True)
lines1,lab1 = ax.get_legend_handles_labels()
lines2,lab2 = ax2.get_legend_handles_labels()
ax.legend(lines1+lines2, lab1+lab2, fontsize=8)

# --- (d) GOR vs SEC_th — 性能权衡 ---
ax=axes[1,1]
sc = ax.scatter(GOR_ok, SEC_th, c=mww_ok, cmap='plasma', s=80, zorder=5)
cb = plt.colorbar(sc, ax=ax); cb.set_label('m_ww [kg/s]')
for i, mw in enumerate(mww_ok):
    if mw % 20 == 0:
        ax.annotate(f'{mw:.0f}', (GOR_ok[i], SEC_th[i]),
                    textcoords="offset points", xytext=(5,3), fontsize=8)
ax.set_xlabel('GOR [-]'); ax.set_ylabel('SEC_th [kJ/kg]')
ax.set_title('(d) GOR ↔ SEC_th 性能权衡\n(颜色=m_ww)', fontsize=11)
ax.grid(True)
# 理论关系：SEC_th = hfg / GOR
gor_range = np.linspace(GOR_ok.min()*0.9, GOR_ok.max()*1.1, 50)
ax.plot(gor_range, 2430/gor_range, 'k--', lw=1.2, label='SEC=2430/GOR 理论线')
ax.legend(fontsize=8)

fig.suptitle('SEC (比能耗) 专项分析 — 单位产水能量消耗', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('result_sec.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ 图2保存：result_sec.png")

# ─────────────────────────────────────────────
#  图 3：收敛性验证
# ─────────────────────────────────────────────
fig, axes = plt.subplots(1,2,figsize=(11,4.5))
fig.patch.set_facecolor('white')

ax=axes[0]
ax.bar(mww_ok, ite_ok, color=C1, alpha=0.8, edgecolor='w')
ax.set_xlabel('m_ww [kg/s]'); ax.set_ylabel('迭代次数')
ax.set_title('SLTD 迭代收敛次数\n（全部工况均已收敛）', fontsize=11); ax.grid(True, axis='y')
for i, v in enumerate(ite_ok):
    ax.text(mww_ok[i], v+5, str(int(v)), ha='center', fontsize=7)

ax=axes[1]
ax.plot(mww_ok, Wsys_ok, '-o', color=C3, **mk)
ax.axhline(0, color='k', lw=1.2, ls='--')
ax.fill_between(mww_ok, Wsys_ok, 0,
                where=(Wsys_ok>0), color='green', alpha=0.25, label='净发电区域')
ax.fill_between(mww_ok, Wsys_ok, 0,
                where=(Wsys_ok<0), color='red', alpha=0.2, label='净耗电区域')
ax.set_xlabel('m_ww [kg/s]'); ax.set_ylabel('W_sys [kW]')
ax.set_title(f'系统净功率 (Hcw={Hcw}m 液压扬程)\nWnet − W_ww泵 − W_cw泵', fontsize=11)
ax.legend(fontsize=9); ax.grid(True)

fig.suptitle('系统收敛性与净功率分析', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('result_convergence.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ 图3保存：result_convergence.png")
