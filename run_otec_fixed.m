% ============================================================
%  run_otec_fixed.m  — 修正 Bug1 + Bug2 后的参数
%    Bug1 修复: mr = 1.5 kg/s （与 UA=500 kW/K 匹配）
%    Bug2 修复: Hcw = 25 m   （实际泵扬程，非管道深度）
%    其余参数与原版完全一致
% ============================================================
clear; clc; close all;
graphics_toolkit('gnuplot');
warning('off','all');

%% 参数
Twwi   = 28;
Tcwi   = 5;
mcw    = 500;
cpw    = 4186;
UAe    = 500e3;    % W/K
UAc    = 500e3;    % W/K
eta_t  = 0.80;
eta_p  = 0.75;
eta_wp = 0.80;
mr     = 1.5;      % ← 修正：原 50 kg/s，现 1.5 kg/s
Hww    = 10;
Hcw    = 25;       % ← 修正：原 800 m，现 25 m（实际扬程）

N_stage = 1;
S0      = 0.035;
NEFe    = 0.02;
NEFc    = 0.02;
hfg_ref = 2.45e6;

mww_vec = (50:50:1000)';
nM      = numel(mww_vec);

fld = {'mww','Teva','Tcon','Twwo','Tcwo', ...
       'h1','h2','h3','h4','s1','s2','s3','s4', ...
       'Wt','Wp','Wnet','Wnet_sys','Qeva','Qcon','eta_cycle', ...
       'Www_pump','Wcw_pump', ...
       'DT','PR','GOR','T_sltd_hot','T_sltd_cold','converged'};
for i = 1:numel(fld), R.(fld{i}) = zeros(nM,1); end

opts = optimset('Display','off','TolFun',1e-10,'TolX',1e-10, ...
                'MaxIter',500,'MaxFunEvals',3000);

x0 = [22; 10; 25; 7];   % 修正初值（避免 Tcon=Tcwo）

fprintf('==================================================\n');
fprintf('  [修正版] 闭式 OTEC + SLTD 耦合系统\n');
fprintf('==================================================\n');
fprintf('  mr=%.1f kg/s  UAe/UAc=%.0f kW/K  Hcw=%.0f m\n\n', mr, UAe/1e3, Hcw);

n_conv = 0;
for k = 1:nM
    mww = mww_vec(k);
    [xk, fval, flag] = fsolve(@(x) residuals_LMTD(x, ...
        Twwi,Tcwi,mww,mcw,mr,UAe,UAc,cpw,eta_t,eta_p), x0, opts);

    if flag <= 0 || norm(fval) > 1e-3
        fprintf('  [未收敛] mww=%4.0f  flag=%d  ||F||=%.2e\n', mww, flag, norm(fval));
        R.converged(k) = 0; continue
    end

    Teva=xk(1); Tcon=xk(2); Twwo=xk(3); Tcwo=xk(4);

    % 增强物理约束（加温度范围限制）
    if Teva>=Twwi || Tcon<=Tcwi || Twwo<=Teva || Tcwo>=Tcon || Teva<0 || Tcon>50
        fprintf('  [越界] mww=%4.0f Teva=%.1f Tcon=%.1f\n', mww, Teva, Tcon);
        R.converged(k) = 0; continue
    end

    [h1,h2,h3,h4,s1,s2,s3,s4,Wt,Wp] = performance_otec(Teva,Tcon,mr,eta_t,eta_p);
    Qeva = mr*(h2-h1)*1e3;
    Qcon = mr*(h3-h4)*1e3;
    Wnet = Wt - Wp;
    eta_cycle = Wnet / Qeva;

    g = 9.81;
    Www_pump = mww*g*Hww / eta_wp;
    Wcw_pump = mcw*g*Hcw / eta_wp;
    Wnet_sys = Wnet - Www_pump - Wcw_pump;

    [~,~,~,~,~,~,~,~,~,~,~,~,T1s,T2s,DT] = ...
        Multistage_in_0423E(N_stage, Twwo, Tcwo, S0, NEFe, NEFc, mww, mcw);
    Q_total = mww*cpw*(Twwo - T1s);   % Bug3 修正: Twwo 而非 Twwi
    if Q_total > 0
        PR  = DT*hfg_ref / Q_total;
        GOR = DT / (Q_total / hfg_ref);
    else
        PR = 0; GOR = 0;
    end

    R.mww(k)=mww;    R.Teva(k)=Teva;  R.Tcon(k)=Tcon;
    R.Twwo(k)=Twwo;  R.Tcwo(k)=Tcwo;
    R.h1(k)=h1; R.h2(k)=h2; R.h3(k)=h3; R.h4(k)=h4;
    R.s1(k)=s1; R.s2(k)=s2; R.s3(k)=s3; R.s4(k)=s4;
    R.Wt(k)=Wt; R.Wp(k)=Wp; R.Wnet(k)=Wnet; R.Wnet_sys(k)=Wnet_sys;
    R.Qeva(k)=Qeva; R.Qcon(k)=Qcon; R.eta_cycle(k)=eta_cycle;
    R.Www_pump(k)=Www_pump; R.Wcw_pump(k)=Wcw_pump;
    R.DT(k)=DT; R.PR(k)=PR; R.GOR(k)=GOR;
    R.T_sltd_hot(k)=T1s; R.T_sltd_cold(k)=T2s;
    R.converged(k) = 1;
    n_conv = n_conv + 1;
    x0 = xk;
end

fprintf('\n扫描完成!  收敛: %d/%d\n\n', n_conv, nM);

%% 结果表格
v = R.converged == 1;
fprintf('======================= [修正版] 工况扫描结果 =======================\n');
fprintf('%5s %6s %6s %6s %6s %7s %9s %6s %8s %7s\n', ...
        'mww','Teva','Tcon','Twwo','Tcwo','Wt(kW)','Wsys(kW)','eta%','DT(kg/s)','GOR');
fprintf('%s\n', repmat('-',1,74));
idx = find(v(:)');
for k = idx
    fprintf('%5.0f %6.2f %6.2f %6.2f %6.2f %7.3f %9.3f %6.3f %8.5f %7.4f\n', ...
        R.mww(k), R.Teva(k), R.Tcon(k), R.Twwo(k), R.Tcwo(k), ...
        R.Wt(k)/1e3, R.Wnet_sys(k)/1e3, R.eta_cycle(k)*100, ...
        R.DT(k), R.GOR(k));
end
fprintf('%s\n\n', repmat('=',1,74));

if any(v)
    k1=idx(1); k_last=idx(end);
    fprintf('--- 关键诊断 ---\n');
    fprintf('  Q_eva      = %.2f kW  (UA×LMTD 极限 ≈ %.0f kW)\n', ...
            R.Qeva(k1)/1e3, UAe*(Twwi-Tcwi)/1e3);
    fprintf('  Wt(透平)   = %.2f kW\n', R.Wt(k1)/1e3);
    fprintf('  Wcw_pump   = %.2f kW  (Hcw=25m)\n', R.Wcw_pump(k1)/1e3);
    fprintf('  Wnet_sys   = %.2f kW\n', R.Wnet_sys(k1)/1e3);
    fprintf('  循环热效率 = %.3f%%\n\n', R.eta_cycle(k1)*100);
end

%% 绘图
if ~any(v)
    fprintf('无收敛工况，跳过绘图。\n');
else
    mwwv = R.mww(v);
    C1=[0.00 0.45 0.74]; C2=[0.85 0.33 0.10];
    C3=[0.47 0.67 0.19]; C4=[0.49 0.18 0.56]; C5=[0.93 0.69 0.13];

    fig1 = figure('Position',[0 0 1280 800]);
    subplot(2,3,1);
    p1=plot(mwwv,R.Wt(v)/1e3,'-o','Color',C1,'LineWidth',1.8,'MarkerSize',5); hold on;
    p2=plot(mwwv,R.Wnet(v)/1e3,'-s','Color',C2,'LineWidth',1.8,'MarkerSize',5);
    p3=plot(mwwv,R.Wnet_sys(v)/1e3,'-d','Color',C3,'LineWidth',1.8,'MarkerSize',5);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('[kW]');
    legend([p1 p2 p3],{'Wt','Wnet','Wsys'},'Location','best');
    title('(a) 发电功率');

    subplot(2,3,2);
    plot(mwwv,R.DT(v),'-o','Color',C2,'LineWidth',1.8,'MarkerSize',5);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('DT [kg/s]');
    title('(b) SLTD 产水量');

    subplot(2,3,3);
    plot(mwwv,R.GOR(v),'-^','Color',C4,'LineWidth',1.8,'MarkerSize',5);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('GOR [-]');
    title('(c) GOR');

    subplot(2,3,4);
    plot(mwwv,R.eta_cycle(v)*100,'-o','Color',C1,'LineWidth',1.8,'MarkerSize',5);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('eta_cycle [%]');
    title('(d) Rankine 循环效率');

    subplot(2,3,5);
    plot(mwwv,R.Wt(v)/1e3,'-o','Color',C1,'LineWidth',1.8,'MarkerSize',5); hold on;
    plot(mwwv,R.Wp(v)/1e3,'--s','Color',C4,'LineWidth',1.5,'MarkerSize',5);
    plot(mwwv,R.Www_pump(v)/1e3,'-.d','Color',C2,'LineWidth',1.5,'MarkerSize',5);
    plot(mwwv,R.Wcw_pump(v)/1e3,'-.^','Color',C5,'LineWidth',1.5,'MarkerSize',5);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('[kW]');
    legend({'Wt','Wp','Www','Wcw'},'Location','best');
    title('(e) 功率分解');

    subplot(2,3,6);
    plot(mwwv,R.Teva(v),'-o','Color',C2,'LineWidth',1.5,'MarkerSize',4); hold on;
    plot(mwwv,R.Tcon(v),'-s','Color',C1,'LineWidth',1.5,'MarkerSize',4);
    plot(mwwv,R.Twwo(v),'--d','Color',C2,'LineWidth',1.3);
    plot(mwwv,R.Tcwo(v),'--^','Color',C1,'LineWidth',1.3);
    plot(mwwv,R.T_sltd_hot(v),':','Color',C4,'LineWidth',1.6);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('[C]');
    legend({'Teva','Tcon','Twwo','Tcwo','T_sltd'},'Location','best');
    title('(f) 关键温度');

    print(fig1,'fig1_fixed.png','-dpng','-r120');
    fprintf('图1已保存: fig1_fixed.png\n');

    fig2 = figure('Position',[0 0 1000 420]);
    subplot(1,2,1);
    plot(mwwv,R.Wnet_sys(v)/1e3,'-o','Color',C1,'LineWidth',2,'MarkerSize',6);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('Wsys [kW]');
    title('系统净功率');
    subplot(1,2,2);
    plot(mwwv,R.DT(v),'-s','Color',C2,'LineWidth',2,'MarkerSize',6);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('DT [kg/s]');
    title('淡水产量');
    print(fig2,'fig2_fixed.png','-dpng','-r120');
    fprintf('图2已保存: fig2_fixed.png\n');
end

save('OTEC_SLTD_fixed.mat','R');
fprintf('结果已保存至 OTEC_SLTD_fixed.mat\n');
