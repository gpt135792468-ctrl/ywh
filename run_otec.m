% ============================================================
%  run_otec.m  — Octave 主运行脚本（与原 easytry.m 参数完全一致）
%  只改了两处语法：optimoptions→optimset，yyaxis→双子图
% ============================================================
clear; clc; close all;
graphics_toolkit('gnuplot');
warning('off', 'all');

%% 边界条件（原版参数，未修改）
Twwi   = 28;
Tcwi   = 5;
mcw    = 500;
cpw    = 4186;
UAe    = 500e3;
UAc    = 500e3;
eta_t  = 0.80;
eta_p  = 0.75;
eta_wp = 0.80;
mr     = 50;
Hww    = 10;
Hcw    = 800;

N_stage = 1;
S0      = 0.035;
NEFe    = 0.02;
NEFc    = 0.02;
hfg_ref = 2.45e6;

%% 扫描热侧流量
mww_vec = (50:50:1000)';
nM      = numel(mww_vec);

fld = {'mww','Teva','Tcon','Twwo','Tcwo', ...
       'h1','h2','h3','h4','s1','s2','s3','s4', ...
       'Wt','Wp','Wnet','Wnet_sys','Qeva','Qcon','eta_cycle', ...
       'Www_pump','Wcw_pump', ...
       'DT','PR','GOR','T_sltd_hot','T_sltd_cold','converged', ...
       'fsolve_flag','fsolve_norm'};
for i = 1:numel(fld), R.(fld{i}) = zeros(nM,1); end

% Octave 用 optimset
opts = optimset('Display','off','TolFun',1e-10,'TolX',1e-10, ...
                'MaxIter',500,'MaxFunEvals',3000);

x0 = [22; 8; 25; 8];

fprintf('==================================================\n');
fprintf('  闭式 OTEC + SLTD 耦合系统 - 热侧流量扫描\n');
fprintf('==================================================\n');
fprintf('参数: mr=%.0f kg/s  UAe/UAc=%.0f kW/K  Hcw=%.0f m\n\n', ...
        mr, UAe/1e3, Hcw);
fprintf('开始扫描 %d 个工况...\n\n', nM);

n_conv = 0;
for k = 1:nM
    mww = mww_vec(k);

    [xk, fval, flag] = fsolve(@(x) residuals_LMTD(x, ...
        Twwi,Tcwi,mww,mcw,mr,UAe,UAc,cpw,eta_t,eta_p), x0, opts);

    R.fsolve_flag(k) = flag;
    R.fsolve_norm(k) = norm(fval);

    if flag <= 0 || norm(fval) > 1e-3
        fprintf('  [未收敛] mww=%4.0f  flag=%d  ||F||=%.2e\n', mww, flag, norm(fval));
        R.converged(k) = 0;
        continue
    end

    Teva = xk(1);  Tcon = xk(2);
    Twwo = xk(3);  Tcwo = xk(4);

    if Teva >= Twwi || Tcon <= Tcwi || Twwo <= Teva || Tcwo >= Tcon
        fprintf('  [违反物理约束] mww=%4.0f  Teva=%.2f Tcon=%.2f Twwo=%.2f Tcwo=%.2f\n', ...
                mww, Teva, Tcon, Twwo, Tcwo);
        R.converged(k) = 0;
        continue
    end

    [h1,h2,h3,h4,s1,s2,s3,s4,Wt,Wp] = performance_otec(Teva, Tcon, mr, eta_t, eta_p);

    Qeva = mr*(h2 - h1)*1e3;
    Qcon = mr*(h3 - h4)*1e3;
    Wnet = Wt - Wp;
    eta_cycle = Wnet / Qeva;

    g = 9.81;
    Www_pump = mww*g*Hww / eta_wp;
    Wcw_pump = mcw*g*Hcw / eta_wp;
    Wnet_sys = Wnet - Www_pump - Wcw_pump;

    [~,~,~,~,~,~,~,~,~,~,~,~,T1s,T2s,DT] = ...
        Multistage_in_0423E(N_stage, Twwo, Tcwo, S0, NEFe, NEFc, mww, mcw);
    Q_total = mww*cpw*(Twwi - T1s);
    if Q_total > 0
        PR  = DT*hfg_ref / Q_total;
        GOR = DT / (Q_total / hfg_ref);
    else
        PR = 0; GOR = 0;
    end

    R.mww(k)=mww;    R.Teva(k)=Teva;  R.Tcon(k)=Tcon;
    R.Twwo(k)=Twwo;  R.Tcwo(k)=Tcwo;
    R.h1(k)=h1;  R.h2(k)=h2;  R.h3(k)=h3;  R.h4(k)=h4;
    R.s1(k)=s1;  R.s2(k)=s2;  R.s3(k)=s3;  R.s4(k)=s4;
    R.Wt(k)=Wt;  R.Wp(k)=Wp;  R.Wnet(k)=Wnet;  R.Wnet_sys(k)=Wnet_sys;
    R.Qeva(k)=Qeva;  R.Qcon(k)=Qcon;  R.eta_cycle(k)=eta_cycle;
    R.Www_pump(k)=Www_pump;  R.Wcw_pump(k)=Wcw_pump;
    R.DT(k)=DT;  R.PR(k)=PR;  R.GOR(k)=GOR;
    R.T_sltd_hot(k)=T1s;  R.T_sltd_cold(k)=T2s;
    R.converged(k) = 1;
    n_conv = n_conv + 1;

    x0 = xk;
end

fprintf('\n扫描完成!  收敛: %d/%d 个工况\n\n', n_conv, nM);

%% 打印结果表格
v = R.converged == 1;
fprintf('========================= 工况扫描结果 =========================\n');
fprintf('%5s %6s %6s %6s %6s %8s %10s %6s %8s %7s %7s\n', ...
        'mww','Teva','Tcon','Twwo','Tcwo','Wnet(kW)','Wsys(kW)','eta%','DT(kg/s)','PR','GOR');
fprintf('%s\n', repmat('-',1,78));
idx = find(v(:)');
for k = idx
    fprintf('%5.0f %6.2f %6.2f %6.2f %6.2f %8.2f %10.2f %6.3f %8.4f %7.4f %7.4f\n', ...
        R.mww(k), R.Teva(k), R.Tcon(k), R.Twwo(k), R.Tcwo(k), ...
        R.Wnet(k)/1e3, R.Wnet_sys(k)/1e3, R.eta_cycle(k)*100, ...
        R.DT(k), R.PR(k), R.GOR(k));
end
fprintf('%s\n\n', repmat('=',1,78));

%% 打印关键数值诊断（帮助发现参数问题）
if any(v)
    k1 = idx(1);
    fprintf('--- 诊断信息 (第一个收敛工况 mww=%.0f kg/s) ---\n', R.mww(k1));
    fprintf('  Q_eva(MJ) = %.2f   [mr*(h2-h1)*1e3]\n', R.Qeva(k1)/1e6);
    fprintf('  UA*LMTD   = %.2f MJ  [理论最大 ≈ UAe × (Twwi-Tcwi) = %.2f MJ]\n', ...
            R.Qeva(k1)/1e6, UAe*(Twwi-Tcwi)/1e6);
    fprintf('  Wt(kW)    = %.2f\n', R.Wt(k1)/1e3);
    fprintf('  Wcw_pump(kW) = %.2f  [= mcw*g*Hcw/eta_wp]\n', R.Wcw_pump(k1)/1e3);
    fprintf('  Wnet_sys(kW) = %.2f\n\n', R.Wnet_sys(k1)/1e3);
end

%% 绘图（只在有收敛结果时绘制）
if ~any(v)
    fprintf('无收敛工况，跳过绘图。\n');
else
    mwwv = R.mww(v);
    C1=[0.00 0.45 0.74]; C2=[0.85 0.33 0.10];
    C3=[0.47 0.67 0.19]; C4=[0.49 0.18 0.56];
    C5=[0.93 0.69 0.13];

    fig1 = figure('Position',[0 0 1280 800]);
    subplot(2,3,1);
    plot(mwwv, R.Wt(v)/1e3,'-o','Color',C1,'LineWidth',1.8,'MarkerSize',5); hold on;
    plot(mwwv, R.Wnet(v)/1e3,'-s','Color',C2,'LineWidth',1.8,'MarkerSize',5);
    plot(mwwv, R.Wnet_sys(v)/1e3,'-d','Color',C3,'LineWidth',1.8,'MarkerSize',5);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('[kW]');
    legend({'Wt','Wnet','Wsys'},'Location','best'); title('(a) 发电功率');

    subplot(2,3,2);
    plot(mwwv, R.DT(v),'-o','Color',C2,'LineWidth',1.8,'MarkerSize',5);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('DT [kg/s]');
    title('(b) SLTD 产水量');

    subplot(2,3,3);
    plot(mwwv, R.GOR(v),'-^','Color',C4,'LineWidth',1.8,'MarkerSize',5);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('GOR [-]');
    title('(c) GOR');

    subplot(2,3,4);
    plot(mwwv, R.PR(v),'-v','Color',C5,'LineWidth',1.8,'MarkerSize',5);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('PR [-]');
    title('(d) PR');

    subplot(2,3,5);
    plot(mwwv, R.eta_cycle(v)*100,'-o','Color',C1,'LineWidth',1.8,'MarkerSize',5);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('eta [%]');
    title('(e) Rankine 循环效率');

    subplot(2,3,6);
    plot(mwwv, R.Teva(v),'-o','Color',C2,'LineWidth',1.5,'MarkerSize',4); hold on;
    plot(mwwv, R.Tcon(v),'-s','Color',C1,'LineWidth',1.5,'MarkerSize',4);
    plot(mwwv, R.Twwo(v),'--d','Color',C2,'LineWidth',1.3);
    plot(mwwv, R.Tcwo(v),'--^','Color',C1,'LineWidth',1.3);
    plot(mwwv, R.T_sltd_hot(v),':','Color',C4,'LineWidth',1.6);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('[C]');
    legend({'Teva','Tcon','Twwo','Tcwo','T_sltd'},'Location','best');
    title('(f) 关键温度');

    print(fig1,'fig1_main.png','-dpng','-r120');
    fprintf('图1已保存: fig1_main.png\n');

    fig2 = figure('Position',[0 0 1000 420]);
    subplot(1,2,1);
    plot(mwwv, R.Wnet_sys(v)/1e3,'-o','Color',C1,'LineWidth',2,'MarkerSize',6);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('Wsys [kW]');
    title('系统净功率');
    subplot(1,2,2);
    plot(mwwv, R.DT(v),'-s','Color',C2,'LineWidth',2,'MarkerSize',6);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('DT [kg/s]');
    title('淡水产量');
    print(fig2,'fig2_combined.png','-dpng','-r120');
    fprintf('图2已保存: fig2_combined.png\n');

    fig3 = figure('Position',[0 0 800 500]);
    plot(mwwv, R.Wp(v)/1e3,'-o','Color',C4,'LineWidth',1.8,'MarkerSize',5); hold on;
    plot(mwwv, R.Www_pump(v)/1e3,'-s','Color',C2,'LineWidth',1.8,'MarkerSize',5);
    plot(mwwv, R.Wcw_pump(v)/1e3,'-d','Color',C1,'LineWidth',1.8,'MarkerSize',5);
    plot(mwwv, (R.Wp(v)+R.Www_pump(v)+R.Wcw_pump(v))/1e3,'--k','LineWidth',2);
    grid on; box on; xlabel('mww [kg/s]'); ylabel('[kW]');
    legend({'Wp(工质泵)','Www_pump','Wcw_pump','总寄生'},'Location','best');
    title('寄生功率分解');
    print(fig3,'fig3_parasitic.png','-dpng','-r120');
    fprintf('图3已保存: fig3_parasitic.png\n');
end

save('OTEC_SLTD_results.mat','R');
fprintf('\n结果已保存至 OTEC_SLTD_results.mat\n');
