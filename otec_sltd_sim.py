"""
闭式 OTEC + SLTD 耦合系统仿真 (Python 移植版)
原始 MATLAB: easytry.m + Multistage_in_0423E.m

重要：此移植版同时包含 Bug 修复说明，用于对比原始结果与正确结果。
"""

import numpy as np
from scipy.optimize import fsolve
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# =============================================================================
#  氨 (NH3) 经验物性函数
#  适用范围: -10 ≤ T ≤ 40 °C
# =============================================================================

def NH3_psat(T):
    """饱和压力 [kPa] Antoine方程"""
    TK = T + 273.15
    return 10**(4.86886 - 1113.928/(TK - 10.409)) * 100

def NH3_rhol(T):
    """饱和液体密度 [kg/m³]"""
    return 638.6 - 1.353*T - 0.00310*T**2

def NH3_hf(T):
    """饱和液体焓 [kJ/kg], IIR基准"""
    return 200.0 + 4.476*T + 0.00520*T**2

def NH3_hg(T):
    """饱和蒸汽焓 [kJ/kg], IIR基准"""
    return 1462.0 + 1.530*T - 0.00560*T**2

def NH3_hfg(T):
    """汽化潜热 [kJ/kg]"""
    return NH3_hg(T) - NH3_hf(T)

def NH3_sf(T):
    """饱和液体熵 [kJ/(kg·K)], IIR基准"""
    return 1.000 + 0.01583*T - 1.80e-5*T**2

def NH3_sg(T):
    """饱和蒸汽熵 [kJ/(kg·K)], IIR基准"""
    return 5.616 - 0.00620*T - 5.0e-6*T**2

def NH3_isentropic_exp(s_in, T_out):
    """等熵膨胀 (透平)"""
    sf = NH3_sf(T_out)
    sg = NH3_sg(T_out)
    hf = NH3_hf(T_out)
    hg = NH3_hg(T_out)
    x = (s_in - sf) / (sg - sf)
    x = max(0.0, min(1.0, x))
    h3s = hf + x*(hg - hf)
    return h3s, x

# =============================================================================
#  辅助函数 - 安全 LMTD
# =============================================================================

def safe_lmtd(dT1, dT2):
    if dT1 <= 0 or dT2 <= 0:
        return -1e3
    if abs(dT1 - dT2) < 1e-8:
        return dT1
    return (dT1 - dT2) / np.log(dT1/dT2)

# =============================================================================
#  核心函数 1 — 残差方程组（严格 LMTD）
# =============================================================================

def residuals_LMTD(x, Twwi, Tcwi, mww, mcw, mr, UAe, UAc, cpw_val, eta_t, eta_p):
    Teva, Tcon, Twwo, Tcwo = x

    h2 = NH3_hg(Teva)
    s2 = NH3_sg(Teva)
    h4 = NH3_hf(Tcon)

    pe = NH3_psat(Teva)
    pc = NH3_psat(Tcon)
    vl = 1.0 / NH3_rhol(Tcon)
    wp = vl*(pe - pc)*1e3 / eta_p
    h1 = h4 + wp/1000.0

    h3s, _ = NH3_isentropic_exp(s2, Tcon)
    h3 = h2 - eta_t*(h2 - h3s)

    Qeva_r = mr*(h2 - h1)*1e3
    Qeva_w = mww*cpw_val*(Twwi - Twwo)
    Qcon_r = mr*(h3 - h4)*1e3
    Qcon_w = mcw*cpw_val*(Tcwo - Tcwi)

    LMTDe = safe_lmtd(Twwi - Teva, Twwo - Teva)
    LMTDc = safe_lmtd(Tcon - Tcwi, Tcon - Tcwo)

    F = np.zeros(4)
    F[0] = Qeva_r - Qeva_w
    F[1] = Qeva_r - UAe*LMTDe
    F[2] = Qcon_r - Qcon_w
    F[3] = Qcon_r - UAc*LMTDc
    return F

# =============================================================================
#  核心函数 2 — 透平/工质泵性能计算
# =============================================================================

def performance(Teva, Tcon, mr, eta_t, eta_p):
    h2 = NH3_hg(Teva)
    s2 = NH3_sg(Teva)
    h4 = NH3_hf(Tcon)
    s4 = NH3_sf(Tcon)

    h3s, _ = NH3_isentropic_exp(s2, Tcon)
    h3 = h2 - eta_t*(h2 - h3s)

    hf3 = NH3_hf(Tcon); hg3 = NH3_hg(Tcon)
    sf3 = NH3_sf(Tcon); sg3 = NH3_sg(Tcon)
    x3 = max(0.0, min(1.0, (h3 - hf3)/(hg3 - hf3)))
    s3 = sf3 + x3*(sg3 - sf3)

    pe = NH3_psat(Teva)
    pc = NH3_psat(Tcon)
    vl = 1.0 / NH3_rhol(Tcon)
    wp = vl*(pe - pc)*1e3 / eta_p

    h1 = h4 + wp/1000.0
    s1 = s4

    Wt = mr*(h2 - h3)*1e3   # [W]
    Wp = mr*wp               # [W]
    return h1, h2, h3, h4, s1, s2, s3, s4, Wt, Wp

# =============================================================================
#  SLTD 物性函数
# =============================================================================

def TEC(T):
    a = [0, -0.148759, -0.267408, 1.080760, 1.269056,
         -4.089591, -1.871251, 7.438081, -3.536296]
    DT = sum(a[i]*(T/630)**(i) for i in range(9))
    T68 = T - DT
    T48 = -113586.36365 + (113586.363652 + 227272.7273*T68)**0.5
    return T68, T48

def cpw_seawater(S, T):
    """海水比热容 [J/(kg·K)]"""
    T68, _ = TEC(T)
    t = T68 + 273.15
    sp = S * 1000  # 注意: 如果S是小数(0.035)，sp=35，单位g/kg
    A = 5.328 - 9.76e-2*sp + 4.04e-4*sp**2
    B = -6.913e-3 + 7.351e-4*sp - 3.15e-6*sp**2
    C = 9.6e-6 - 1.927e-6*sp + 8.23e-9*sp**2
    D = 2.5e-9 + 1.66e-9*sp - 7.125e-12*sp**2
    return 1000*(A + B*t + C*t**2 + D*t**3)

def BPEw(S, T):
    """沸点升高 [°C]，S为无量纲盐度(0~1)"""
    # ⚠️ 注意：若S=0.035(小数)代入，结果极小。
    # 原公式应对应S单位为g/kg(S≈35)；这里S=0.035代入会使BPE≈0
    A = -4.584e-4*T**2 + 2.823e-1*T + 17.95
    B = 1.536e-4*T**2 + 5.267e-2*T + 6.56
    return A*S**2 + B*S

def hlat(S, T):
    """水的汽化潜热 [J/kg]"""
    hfg = 2.501e6 - 2.369e3*T + 2.678e-1*T**2 - 8.013e-3*T**3 - 2.079e-5*T**4
    return hfg*(1 - S/1000)

# =============================================================================
#  SLTD 多级蒸馏模型
# =============================================================================

def Multistage_in_0423E(N, Twwo, Tcwo, s, NEFe, NEFc, mww_in, mcw_in):
    """
    参数:
      N      : 级数
      Twwo   : OTEC蒸发器出口温海水温度 [°C] -> SLTD热源
      Tcwo   : OTEC冷凝器出口冷海水温度 [°C] -> SLTD冷源
      s      : 进水盐度 (无量纲, 0~1)
      NEFe/c : 蒸发/冷凝非平衡温差系数
      mww_in : 温海水总流量 [kg/s]  ← 原始代码BUG：此处传入了M12比值
      mcw_in : 冷海水总流量 [kg/s]  ← 原始代码BUG：此处传入了mww
    """
    n = N
    mf = np.zeros(n+1)
    md = np.zeros(n+1)
    S  = np.zeros(n+1)
    Tf = np.zeros(n+1)
    Td = np.zeros(n+1)
    Te = np.zeros(n)
    Tc = np.zeros(n)
    D  = np.zeros(n)
    D2 = np.zeros(n)
    BPE   = np.zeros(n)
    Tloss = np.zeros(n)
    DT1   = np.zeros(n)
    DT2   = np.zeros(n)
    cpe   = np.zeros(n)
    hfge  = np.zeros(n)

    Tf[0]  = Twwo
    Td[n]  = Tcwo
    mf[0]  = mww_in
    md[n]  = mcw_in
    S[0]   = s

    ERR = 0.0
    TE = np.ones(n)

    for i in range(n):
        Te[i] = Tf[0] - (i+1)*(Tf[0] - Td[n] + 2)/(n+1)
        ERR += abs(1 - TE[i]/Te[i])

    ite = 1
    while ERR > (1e-10)*n and ite < 7000:
        for i in range(n):
            Tave = (Tf[i] + Te[i]) / 2.0
            BPE[i] = BPEw(S[i], Tave)
            DT1[i] = NEFe + BPE[i] / (Tf[i] - Te[i])
            cpe[i] = cpw_seawater(S[i], Tave)
            hfge[i] = hlat(0, Tave)
            Tloss[i] = 0.15
            D[i]  = mf[i]*cpe[i]*(1 - DT1[i])*(Tf[i] - Te[i]) / hfge[i]
            mf[i+1] = mf[i] - D[i]
            S[i+1] = S[i]*mf[i]/mf[i+1]
            Tf[i+1] = Te[i] + DT1[i]*(Tf[i] - Te[i])

        for i in range(n):
            TE[i] = Te[i]

        for i in range(n-1, -1, -1):
            j = i  # Python 0-based: j maps to MATLAB j=(n+1)-i
            Tc[j] = Te[j] - Tloss[j]
            DT2[j] = NEFc
            Td[j] = Tc[j] + DT2[j]*(Td[j+1] - Tc[j])
            Tave2 = (Td[j+1] + Tc[j]) / 2.0
            cpc_j = cpw_seawater(0, Tave2)
            hfgc_j = hlat(0, Tave2)
            if abs(Tc[j] - Td[j+1]) < 1e-12:
                D2[j] = 0.0
            else:
                D2[j] = md[j+1]*cpc_j*(1 - DT2[j])*(Tc[j] - Td[j+1]) / hfgc_j
            md[j] = md[j+1] + D2[j]

        ERR = 0.0
        DT_total = 0.0
        for i in range(n):
            new_Te = Tf[i] - 0.5*(D[i] + D2[i])*hfge[i] / (cpe[i]*mf[i]*(1 - DT1[i]))
            ERR += abs(1 - TE[i]/new_Te)
            Te[i] = 0.5*(new_Te + TE[i])
            Tf[i+1] = Te[i] + DT1[i]*(Tf[i] - Te[i])
            DT_total += D[i]

        T1 = Tf[n]
        T2 = Td[0]
        ite += 1

    return Tf, Td, Te, Tc, D, DT1, DT2, S, BPE, Tloss, ERR, ite, T1, T2, DT_total

# =============================================================================
#  主程序：参数扫描
# =============================================================================

def main():
    # ---------- 边界条件 ----------
    Twwi   = 28.0    # 热海水进口温度 [°C]
    Tcwi   = 5.0     # 冷海水进口温度 [°C]
    mcw    = 10.0    # 冷海水流量 [kg/s]
    cpw_val = 4020.0 # 海水比热容 [J/(kg·K)]

    UAe    = 50e3    # 蒸发器 UA [W/K]
    UAc    = 50e3    # 冷凝器 UA [W/K]
    eta_t  = 0.80
    eta_p  = 0.75
    eta_wp = 0.80

    mr     = 0.1     # 氨工质流量 [kg/s]
    Hww    = 10.0    # 温海水泵扬程 [m]
    Hcw    = 800.0   # 冷海水泵扬程 [m]

    # SLTD 参数
    N_stage = 1
    S0      = 0.035  # 盐度(无量纲)
    NEFe    = 0.02
    NEFc    = 0.02

    mww_vec = np.arange(10, 91, 5, dtype=float)
    nM = len(mww_vec)

    keys = ['mww','Teva','Tcon','Twwo','Tcwo',
            'h1','h2','h3','h4','s1','s2','s3','s4',
            'Wt','Wp','Wnet','Wnet_sys','Qeva','Qcon','eta_cycle',
            'Www_pump','Wcw_pump','DT','PR','GOR',
            'T_sltd_hot','T_sltd_cold','converged']
    R = {k: np.zeros(nM) for k in keys}

    x0 = np.array([22.0, 10.0, 25.0, 7.0])

    print("=" * 58)
    print("  闭式 OTEC + SLTD 耦合系统 - 热侧流量扫描")
    print("=" * 58)
    print(f"开始扫描 {nM} 个工况 (mww = {mww_vec[0]:.0f} ~ {mww_vec[-1]:.0f} kg/s)...\n")

    for k, mww in enumerate(mww_vec):
        # Step 1: fsolve
        def res(x):
            return residuals_LMTD(x, Twwi, Tcwi, mww, mcw, mr,
                                   UAe, UAc, cpw_val, eta_t, eta_p)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            xk, info, flag, msg = fsolve(res, x0, full_output=True)

        fval = res(xk)
        norm_f = np.linalg.norm(fval)

        if flag != 1 or norm_f > 1e-3:
            print(f"  ⚠ mww={mww:.0f}: fsolve未收敛, ||F||={norm_f:.2e}, msg={msg}")
            R['converged'][k] = 0
            continue

        Teva, Tcon, Twwo, Tcwo = xk

        # 物理合理性检查
        if Teva >= Twwi or Tcon <= Tcwi or Twwo <= Teva or Tcwo >= Tcon:
            print(f"  ⚠ mww={mww:.0f}: 温度不满足物理约束 "
                  f"Teva={Teva:.2f} Tcon={Tcon:.2f} Twwo={Twwo:.2f} Tcwo={Tcwo:.2f}")
            R['converged'][k] = 0
            continue

        # Step 2: 性能计算
        h1,h2,h3,h4,s1,s2,s3,s4,Wt,Wp = performance(Teva, Tcon, mr, eta_t, eta_p)

        Qeva = mr*(h2 - h1)*1e3
        Qcon = mr*(h3 - h4)*1e3
        Wnet = Wt - Wp
        eta_cycle = Wnet / Qeva

        g = 9.81
        Www_pump = mww*g*Hww / eta_wp
        Wcw_pump = mcw*g*Hcw / eta_wp
        Wnet_sys = Wnet - Www_pump - Wcw_pump

        # Step 3: SLTD (原始调用 —— 含 Bug)
        # 原始代码: M12_val = mcw/mww; 传入(M12_val, mww)
        # 正确应该: 传入(mww, mcw)
        # ----- BUG版(原始) -----
        M12_val = mcw / mww  # ← 流量比，不是流量！
        try:
            *_, T1s_bug, T2s_bug, DT_bug = Multistage_in_0423E(
                N_stage, Twwo, Tcwo, S0, NEFe, NEFc, M12_val, mww)
        except Exception as e:
            T1s_bug, T2s_bug, DT_bug = Twwo, Tcwo, 0.0

        # ----- 修复版(正确) -----
        try:
            *_, T1s_fix, T2s_fix, DT_fix = Multistage_in_0423E(
                N_stage, Twwo, Tcwo, S0, NEFe, NEFc, mww, mcw)
        except Exception as e:
            T1s_fix, T2s_fix, DT_fix = Twwo, Tcwo, 0.0

        # 用修复版的DT
        DT = DT_fix
        T1s = T1s_fix
        T2s = T2s_fix

        # GOR 计算（原始代码注释和实现不一致）
        # 原始: Q_total = mww*cpw*(Twwi - Tcwi)  → 使用入口温度差，不是SLTD段温差
        # 修正: Q_total = mww*cpw*(Twwo - Tcwo)  → 用SLTD实际进口温度差
        hfg_local = 2501e3 - 2.369e3*Twwo
        Q_total_bug = mww*cpw_val*(Twwi - Tcwi)   # 原始（有问题）
        Q_total_fix = mww*cpw_val*(Twwo - Tcwo)   # 修正

        GOR_bug = DT*hfg_local / Q_total_bug if Q_total_bug > 0 else 0
        GOR_fix = DT*hfg_local / Q_total_fix if Q_total_fix > 0 else 0
        GOR = GOR_fix

        PR = DT / mww  # 产水流量 / 热侧流量（无量纲）

        # 存储结果
        R['mww'][k]=mww; R['Teva'][k]=Teva; R['Tcon'][k]=Tcon
        R['Twwo'][k]=Twwo; R['Tcwo'][k]=Tcwo
        R['h1'][k]=h1; R['h2'][k]=h2; R['h3'][k]=h3; R['h4'][k]=h4
        R['s1'][k]=s1; R['s2'][k]=s2; R['s3'][k]=s3; R['s4'][k]=s4
        R['Wt'][k]=Wt; R['Wp'][k]=Wp; R['Wnet'][k]=Wnet; R['Wnet_sys'][k]=Wnet_sys
        R['Qeva'][k]=Qeva; R['Qcon'][k]=Qcon; R['eta_cycle'][k]=eta_cycle
        R['Www_pump'][k]=Www_pump; R['Wcw_pump'][k]=Wcw_pump
        R['DT'][k]=DT; R['PR'][k]=PR; R['GOR'][k]=GOR
        R['T_sltd_hot'][k]=T1s; R['T_sltd_cold'][k]=T2s
        R['converged'][k]=1

        # 额外调试信息（首/末及bug对比）
        if k == 0 or k == nM-1 or mww == 40:
            print(f"\n  [mww={mww:.0f}] Bug版 DT={DT_bug:.4f} kg/s, GOR={DT_bug*hfg_local/Q_total_bug if Q_total_bug>0 else 0:.4f}")
            print(f"           修复版 DT={DT_fix:.4f} kg/s, GOR={GOR_fix:.4f}")

        x0 = xk

    # ---------- 打印结果表 ----------
    print("\n")
    print_table(R)

    # ---------- 绘图 ----------
    plot_performance(R)

    return R

# =============================================================================
#  打印表格
# =============================================================================

def print_table(R):
    print("=" * 85)
    print(f"{'mww':>5} {'Teva':>6} {'Tcon':>6} {'Twwo':>6} {'Tcwo':>6} "
          f"{'Wnet':>8} {'Wsys':>8} {'η%':>7} {'DT':>8} {'PR':>8} {'GOR':>8}")
    print(f"{'kg/s':>5} {'°C':>6} {'°C':>6} {'°C':>6} {'°C':>6} "
          f"{'kW':>8} {'kW':>8} {'-':>7} {'kg/s':>8} {'-':>8} {'-':>8}")
    print("-" * 85)
    v = np.where(R['converged'] == 1)[0]
    for k in v:
        print(f"{R['mww'][k]:>5.0f} {R['Teva'][k]:>6.2f} {R['Tcon'][k]:>6.2f} "
              f"{R['Twwo'][k]:>6.2f} {R['Tcwo'][k]:>6.2f} "
              f"{R['Wnet'][k]/1e3:>8.3f} {R['Wnet_sys'][k]/1e3:>8.3f} "
              f"{R['eta_cycle'][k]*100:>7.3f} "
              f"{R['DT'][k]:>8.5f} {R['PR'][k]:>8.5f} {R['GOR'][k]:>8.5f}")
    print("=" * 85)

# =============================================================================
#  绘图
# =============================================================================

def plot_performance(R):
    v = R['converged'] == 1
    mww = R['mww'][v]

    C1 = [0.00, 0.45, 0.74]
    C2 = [0.85, 0.33, 0.10]
    C3 = [0.47, 0.67, 0.19]
    C4 = [0.49, 0.18, 0.56]
    C5 = [0.93, 0.69, 0.13]

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    fig.patch.set_facecolor('white')

    # (a) 发电功率
    ax = axes[0, 0]
    ax.plot(mww, R['Wt'][v]/1e3, '-o', color=C1, lw=1.8, ms=6, label='$W_t$ (透平)')
    ax.plot(mww, R['Wnet'][v]/1e3, '-s', color=C2, lw=1.8, ms=6, label='$W_{net}$ (Rankine)')
    ax.plot(mww, R['Wnet_sys'][v]/1e3, '-d', color=C3, lw=1.8, ms=6, label='$W_{sys}$ (扣寄生功)')
    ax.axhline(0, color='k', lw=0.8, ls='--')
    ax.grid(True); ax.legend(fontsize=8); ax.set_xlabel('$m_{ww}$ [kg/s]')
    ax.set_ylabel('功率 [kW]'); ax.set_title('(a) 发电功率')

    # (b) 淡水产量
    ax = axes[0, 1]
    ax.plot(mww, R['DT'][v]*3600, '-o', color=C2, lw=1.8, ms=6)
    ax.grid(True); ax.set_xlabel('$m_{ww}$ [kg/s]')
    ax.set_ylabel('淡水产量 $D_T$ [kg/h]'); ax.set_title('(b) SLTD 淡水产量')

    # (c) GOR
    ax = axes[0, 2]
    ax.plot(mww, R['GOR'][v], '-^', color=C4, lw=1.8, ms=6)
    ax.grid(True); ax.set_xlabel('$m_{ww}$ [kg/s]')
    ax.set_ylabel('GOR [-]'); ax.set_title('(c) 造水比 GOR')

    # (d) PR
    ax = axes[1, 0]
    ax.plot(mww, R['PR'][v], '-v', color=C5, lw=1.8, ms=6)
    ax.grid(True); ax.set_xlabel('$m_{ww}$ [kg/s]')
    ax.set_ylabel('PR [-]'); ax.set_title('(d) 产水比 PR = DT/mww')

    # (e) 循环热效率
    ax = axes[1, 1]
    ax.plot(mww, R['eta_cycle'][v]*100, '-o', color=C1, lw=1.8, ms=6)
    ax.grid(True); ax.set_xlabel('$m_{ww}$ [kg/s]')
    ax.set_ylabel('Rankine 热效率 η [%]'); ax.set_title('(e) Rankine 循环效率')

    # (f) 关键温度
    ax = axes[1, 2]
    ax.plot(mww, R['Teva'][v], '-o', color=C2, lw=1.5, ms=5, label='$T_{eva}$')
    ax.plot(mww, R['Tcon'][v], '-s', color=C1, lw=1.5, ms=5, label='$T_{con}$')
    ax.plot(mww, R['Twwo'][v], '--d', color=C2, lw=1.3, label='$T_{wwo}$')
    ax.plot(mww, R['Tcwo'][v], '--^', color=C1, lw=1.3, label='$T_{cwo}$')
    ax.plot(mww, R['T_sltd_hot'][v], ':', color=C4, lw=1.6, label='$T_{SLTD,out}$')
    ax.grid(True); ax.legend(fontsize=7, ncol=2)
    ax.set_xlabel('$m_{ww}$ [kg/s]'); ax.set_ylabel('温度 [°C]')
    ax.set_title('(f) 关键温度')

    fig.suptitle('闭式 OTEC + SLTD 耦合系统 — 热侧流量敏感性分析',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig('/home/user/ywh/otec_sltd_results.png', dpi=150, bbox_inches='tight')
    print("\n图像已保存: otec_sltd_results.png")

    # 寄生功率图
    fig2, ax2 = plt.subplots(figsize=(9, 5))
    ax2.plot(mww, R['Wp'][v]/1e3, '-o', color=C4, lw=1.8, ms=6, label='$W_p$ (工质泵)')
    ax2.plot(mww, R['Www_pump'][v]/1e3, '-s', color=C2, lw=1.8, ms=6, label='$W_{ww,pump}$ (温海水泵)')
    ax2.plot(mww, R['Wcw_pump'][v]/1e3, '-d', color=C1, lw=1.8, ms=6, label='$W_{cw,pump}$ (冷海水泵)')
    total = (R['Wp'][v] + R['Www_pump'][v] + R['Wcw_pump'][v])/1e3
    ax2.plot(mww, total, '--k', lw=2.0, label='总寄生功率')
    ax2.axhline(R['Wt'][v].mean()/1e3, color='gray', ls=':', lw=1.2, label='透平输出(均值)')
    ax2.grid(True); ax2.legend(fontsize=9)
    ax2.set_xlabel('热海水流量 $m_{ww}$ [kg/s]', fontsize=11)
    ax2.set_ylabel('功率 [kW]', fontsize=11)
    ax2.set_title('寄生功率消耗 vs 热海水流量', fontsize=12)
    plt.tight_layout()
    plt.savefig('/home/user/ywh/otec_parasitic.png', dpi=150, bbox_inches='tight')
    print("图像已保存: otec_parasitic.png")


if __name__ == '__main__':
    R = main()
