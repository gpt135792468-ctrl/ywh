clear; clc; close all;
graphics_toolkit('gnuplot');   % Octave headless 绘图后端

%% 边界条件
%海水参数
Twwi   = 28;        % 热海水进口温度 [°C]
Tcwi   = 5;         % 冷海水进口温度 [°C]
mcw    = 500;       % 冷海水流量 [kg/s]（固定）
cpw    = 4186;      % 海水比热容 [J/(kg·K)]
%OTEC 设备参数
UAe    = 500e3;     % 蒸发器总传热系数×面积 [W/K]
UAc    = 500e3;     % 冷凝器总传热系数×面积 [W/K]
eta_t  = 0.80;      % 透平等熵效率 [-]
eta_p  = 0.75;      % 工质泵效率 [-]
eta_wp = 0.80;      % 海水泵效率 [-]
%工质流量
mr     = 50;       % 氨工质流量 [kg/s]
% 海水泵扬程
Hww    = 10;        % 温海水抽水高度 [m]
Hcw    = 800;       % 冷海水抽水高度 [m]

% --- SLTD 参数 ---
N_stage = 1;
S0      = 0.035;
NEFe    = 0.02;
NEFc    = 0.02;
hfg_ref = 2.45e6;   % 参考潜热 [J/kg]

%% ========== 2. 扫描热侧流量 ==========

mww_vec = (50:50:1000)';
nM      = numel(mww_vec);

% --- 预分配 ---
fld = {'mww','Teva','Tcon','Twwo','Tcwo', ...
       'h1','h2','h3','h4','s1','s2','s3','s4', ...
       'Wt','Wp','Wnet','Wnet_sys','Qeva','Qcon','eta_cycle', ...
       'Www_pump','Wcw_pump', ...
       'DT','PR','GOR','T_sltd_hot','T_sltd_cold','converged'};
for i = 1:numel(fld), R.(fld{i}) = zeros(nM,1); end

% ---- Octave 兼容：optimset 代替 optimoptions ----
opts = optimset('Display','off', ...
       'TolFun',1e-10,'TolX',1e-10, ...
       'MaxIter',500,'MaxFunEvals',3000);

x0 = [22; 8; 25; 8];      % 初值

fprintf('==================================================\n');
fprintf('  闭式 OTEC + SLTD 耦合系统 - 热侧流量扫描\n');
fprintf('==================================================\n');
fprintf('开始扫描 %d 个工况 (mww = %g ~ %g kg/s)...\n\n', ...
    nM, mww_vec(1), mww_vec(end));

for k = 1:nM
    mww = mww_vec(k);

    % ----- Step 1: 求解 4×4 非线性方程组 -----
    [xk, fval, flag] = fsolve(@(x) residuals_LMTD(x, ...
        Twwi,Tcwi,mww,mcw,mr,UAe,UAc,cpw,eta_t,eta_p), x0, opts);

    if flag <= 0 || norm(fval) > 1e-3
        warning('mww = %.0f: fsolve 未收敛, ||F|| = %.2e', mww, norm(fval));
        R.converged(k) = 0;
        continue
    end

    Teva = xk(1);  Tcon = xk(2);
    Twwo = xk(3);  Tcwo = xk(4);

    % 物理合理性检查
    if Teva >= Twwi || Tcon <= Tcwi || Twwo <= Teva || Tcwo >= Tcon
        warning('mww = %.0f: 温度不满足物理约束', mww);
        R.converged(k) = 0;
        continue
    end

    % ----- Step 2: 性能计算 -----
    [h1,h2,h3,h4,s1,s2,s3,s4,Wt,Wp] = ...
        performance(Teva, Tcon, mr, eta_t, eta_p);

    Qeva = mr*(h2 - h1)*1e3;          % [W]
    Qcon = mr*(h3 - h4)*1e3;          % [W]
    Wnet = Wt - Wp;                   % [W]
    eta_cycle = Wnet / Qeva;

    g = 9.81;
    Www_pump = mww*g*Hww / eta_wp;    % [W]
    Wcw_pump = mcw*g*Hcw / eta_wp;    % [W]
    Wnet_sys = Wnet - Www_pump - Wcw_pump;

    % ----- Step 3: SLTD -----
    [~,~,~,~,~,~,~,~,~,~,~,~,T1s,T2s,DT] = ...
    Multistage_in_0423E(N_stage, Twwo, Tcwo, S0, NEFe, NEFc, mww, mcw);
    Q_total = mww*cpw*(Twwi - T1s);
    if Q_total > 0
        PR  = DT*hfg_ref / Q_total;
        GOR = DT / (Q_total / hfg_ref);
    else
        PR = 0; GOR = 0;
    end

    % ----- Step 4: 存储 -----
    R.mww(k)=mww; R.Teva(k)=Teva; R.Tcon(k)=Tcon;
    R.Twwo(k)=Twwo; R.Tcwo(k)=Tcwo;
    R.h1(k)=h1; R.h2(k)=h2; R.h3(k)=h3; R.h4(k)=h4;
    R.s1(k)=s1; R.s2(k)=s2; R.s3(k)=s3; R.s4(k)=s4;
    R.Wt(k)=Wt; R.Wp(k)=Wp; R.Wnet(k)=Wnet; R.Wnet_sys(k)=Wnet_sys;
    R.Qeva(k)=Qeva; R.Qcon(k)=Qcon; R.eta_cycle(k)=eta_cycle;
    R.Www_pump(k)=Www_pump; R.Wcw_pump(k)=Wcw_pump;
    R.DT(k)=DT; R.PR(k)=PR; R.GOR(k)=GOR;
    R.T_sltd_hot(k)=T1s; R.T_sltd_cold(k)=T2s;
    R.converged(k) = 1;

    x0 = xk;   % 用当前解做下一工况的初值

    if mod(k, 5) == 0
        fprintf('  进度: %d/%d (%.0f%%)\n', k, nM, 100*k/nM);
    end
end
fprintf('扫描完成!\n\n');

%% ========== 3. 打印结果表格 ==========
print_table(R);

%% ========== 4. 绘制性能评估图 ==========
plot_performance(R);

%% ========== 5. 保存数据 ==========
save('OTEC_SLTD_results.mat', 'R');
fprintf('结果已保存至 OTEC_SLTD_results.mat\n');


%% =========================================================================
%       核心函数 1 — 残差方程组（严格 LMTD）
% =========================================================================
function F = residuals_LMTD(x, Twwi,Tcwi,mww,mcw,mr,UAe,UAc,cpw,eta_t,eta_p)
    Teva = x(1);  Tcon = x(2);
    Twwo = x(3);  Tcwo = x(4);

    % ---------- 工质侧焓 ----------
    h2 = NH3_hg(Teva);
    s2 = NH3_sg(Teva);
    h4 = NH3_hf(Tcon);

    % 工质泵 (4→1)
    pe = NH3_psat(Teva);
    pc = NH3_psat(Tcon);
    vl = 1 / NH3_rhol(Tcon);
    wp = vl*(pe - pc)*1e3 / eta_p;
    h1 = h4 + wp/1000;
    % 透平 (2→3)
    [h3s, ~] = NH3_isentropic_exp(s2, Tcon);
    h3 = h2 - eta_t*(h2 - h3s);

    % ---------- 热量 ----------
    Qeva_r = mr*(h2 - h1)*1e3;
    Qeva_w = mww*cpw*(Twwi - Twwo);
    Qcon_r = mr*(h3 - h4)*1e3;
    Qcon_w = mcw*cpw*(Tcwo - Tcwi);

    % ---------- LMTD ----------
    LMTDe = safe_lmtd(Twwi - Teva, Twwo - Teva);
    LMTDc = safe_lmtd(Tcon - Tcwi, Tcon - Tcwo);

    F = zeros(4,1);
    F(1) = Qeva_r - Qeva_w;
    F(2) = Qeva_r - UAe*LMTDe;
    F(3) = Qcon_r - Qcon_w;
    F(4) = Qcon_r - UAc*LMTDc;
end

%% =========================================================================
%       核心函数 2 — 透平 / 工质泵 性能计算
% =========================================================================
function [h1,h2,h3,h4,s1,s2,s3,s4,Wt,Wp] = ...
         performance(Teva, Tcon, mr, eta_t, eta_p)

    h2 = NH3_hg(Teva);   s2 = NH3_sg(Teva);
    h4 = NH3_hf(Tcon);   s4 = NH3_sf(Tcon);

    [h3s, ~] = NH3_isentropic_exp(s2, Tcon);
    h3 = h2 - eta_t*(h2 - h3s);

    hf3 = NH3_hf(Tcon); hg3 = NH3_hg(Tcon);
    sf3 = NH3_sf(Tcon); sg3 = NH3_sg(Tcon);
    x3  = max(0, min(1, (h3-hf3)/(hg3-hf3)));
    s3  = sf3 + x3*(sg3 - sf3);

    pe = NH3_psat(Teva);  pc = NH3_psat(Tcon);
    vl = 1/NH3_rhol(Tcon);
    wp = vl*(pe - pc)*1e3 / eta_p;

    h1 = h4 + wp/1000;
    s1 = s4;

    Wt = mr*(h2 - h3)*1e3;
    Wp = mr*wp;
end

%% =========================================================================
%       打印表格
% =========================================================================
function print_table(R)
    fprintf('====================== 工况扫描结果 ======================\n');
    fprintf('%5s %6s %6s %6s %6s %8s %8s %7s %7s %7s %7s\n', ...
        'mww','Teva','Tcon','Twwo','Tcwo','Wnet','Wsys','eta%','DT','PR','GOR');
    fprintf('%5s %6s %6s %6s %6s %8s %8s %7s %7s %7s %7s\n', ...
        'kg/s','C','C','C','C','kW','kW','-','kg/s','-','-');
    fprintf('------------------------------------------------------------------------\n');
    v = R.converged == 1;
    idx = find(v(:)');
    for k = idx
        fprintf('%5.0f %6.2f %6.2f %6.2f %6.2f %8.2f %8.2f %7.3f %7.4f %7.4f %7.4f\n', ...
            R.mww(k),R.Teva(k),R.Tcon(k),R.Twwo(k),R.Tcwo(k), ...
            R.Wnet(k)/1e3, R.Wnet_sys(k)/1e3, R.eta_cycle(k)*100, ...
            R.DT(k), R.PR(k), R.GOR(k));
    end
    fprintf('========================================================================\n\n');
end

%% =========================================================================
%       绘图函数 — 性能评估图
% =========================================================================
function plot_performance(R)
    v = R.converged == 1;
    mww = R.mww(v);

    C1 = [0.00 0.45 0.74];
    C2 = [0.85 0.33 0.10];
    C3 = [0.47 0.67 0.19];
    C4 = [0.49 0.18 0.56];
    C5 = [0.93 0.69 0.13];

    % ============= 主性能图 =============
    fig1 = figure('Position',[60 60 1280 800]);

    subplot(2,3,1);
    p1 = plot(mww, R.Wt(v)/1e3, '-o','Color',C1,'LineWidth',1.8,'MarkerSize',6); hold on;
    p2 = plot(mww, R.Wnet(v)/1e3, '-s','Color',C2,'LineWidth',1.8,'MarkerSize',6);
    p3 = plot(mww, R.Wnet_sys(v)/1e3, '-d','Color',C3,'LineWidth',1.8,'MarkerSize',6);
    grid on; box on;
    xlabel('热海水流量 mww [kg/s]'); ylabel('功率 [kW]');
    legend([p1 p2 p3],{'Wt (透平)','Wnet (Rankine)','Wsys (扣寄生功)'},'Location','best');
    title('(a) 发电功率 vs mww');

    subplot(2,3,2);
    plot(mww, R.DT(v), '-o','Color',C2,'LineWidth',1.8,'MarkerSize',6);
    grid on; box on;
    xlabel('热海水流量 mww [kg/s]'); ylabel('总产水量 DT [kg/s]');
    title('(b) SLTD 淡水产量 vs mww');

    subplot(2,3,3);
    plot(mww, R.GOR(v), '-^','Color',C4,'LineWidth',1.8,'MarkerSize',6);
    grid on; box on;
    xlabel('热海水流量 mww [kg/s]'); ylabel('GOR [-]');
    title('(c) 造水比 GOR vs mww');

    subplot(2,3,4);
    plot(mww, R.PR(v), '-v','Color',C5,'LineWidth',1.8,'MarkerSize',6);
    grid on; box on;
    xlabel('热海水流量 mww [kg/s]'); ylabel('PR [-]');
    title('(d) 性能比 PR vs mww');

    subplot(2,3,5);
    plot(mww, R.eta_cycle(v)*100, '-o','Color',C1,'LineWidth',1.8,'MarkerSize',6);
    grid on; box on;
    xlabel('热海水流量 mww [kg/s]'); ylabel('Rankine 循环热效率 [%]');
    title('(e) Rankine 循环效率 vs mww');

    subplot(2,3,6);
    plot(mww, R.Teva(v), '-o','Color',C2,'LineWidth',1.5,'MarkerSize',5); hold on;
    plot(mww, R.Tcon(v), '-s','Color',C1,'LineWidth',1.5,'MarkerSize',5);
    plot(mww, R.Twwo(v), '--d','Color',C2,'LineWidth',1.3);
    plot(mww, R.Tcwo(v), '--^','Color',C1,'LineWidth',1.3);
    plot(mww, R.T_sltd_hot(v), ':','Color',C4,'LineWidth',1.6);
    grid on; box on;
    xlabel('热海水流量 mww [kg/s]'); ylabel('温度 [C]');
    legend({'Teva','Tcon','Twwo','Tcwo','T_SLTD_out'},'Location','best');
    title('(f) 关键温度 vs mww');

    print(fig1, 'fig1_main_performance.png', '-dpng', '-r150');
    fprintf('图1已保存: fig1_main_performance.png\n');

    % ============= GOR & PR 详细图 =============
    fig2 = figure('Position',[120 120 1100 450]);
    subplot(1,2,1);
    plot(mww, R.GOR(v), '-^','Color',C4,'LineWidth',2.0,'MarkerSize',7);
    grid on; box on;
    xlabel('热海水流量 mww [kg/s]'); ylabel('GOR [-]');
    title('造水比 GOR vs 热海水流量');

    subplot(1,2,2);
    plot(mww, R.PR(v), '-v','Color',C5,'LineWidth',2.0,'MarkerSize',7);
    grid on; box on;
    xlabel('热海水流量 mww [kg/s]'); ylabel('PR [-]');
    title('性能比 PR vs 热海水流量');

    print(fig2, 'fig2_GOR_PR.png', '-dpng', '-r150');
    fprintf('图2已保存: fig2_GOR_PR.png\n');

    % ============= 发电功率详细图 =============
    fig3 = figure('Position',[180 180 800 500]);
    plot(mww, R.Wt(v)/1e3,       '-o','Color',C1,'LineWidth',2.0,'MarkerSize',7); hold on;
    plot(mww, R.Wnet(v)/1e3,     '-s','Color',C2,'LineWidth',2.0,'MarkerSize',7);
    plot(mww, R.Wnet_sys(v)/1e3, '-d','Color',C3,'LineWidth',2.0,'MarkerSize',7);
    grid on; box on;
    xlabel('热海水流量 mww [kg/s]'); ylabel('功率 [kW]');
    legend({'Wt (透平输出)','Wnet (Rankine净功)','Wsys (系统净功)'},'Location','best');
    title('不同热海水流量下的发电功率');

    print(fig3, 'fig3_power.png', '-dpng', '-r150');
    fprintf('图3已保存: fig3_power.png\n');

    % ============= 寄生功率分解 =============
    fig4 = figure('Position',[240 240 800 500]);
    plot(mww, R.Wp(v)/1e3,       '-o','Color',C4,'LineWidth',1.8,'MarkerSize',6); hold on;
    plot(mww, R.Www_pump(v)/1e3, '-s','Color',C2,'LineWidth',1.8,'MarkerSize',6);
    plot(mww, R.Wcw_pump(v)/1e3, '-d','Color',C1,'LineWidth',1.8,'MarkerSize',6);
    plot(mww, (R.Wp(v)+R.Www_pump(v)+R.Wcw_pump(v))/1e3, '--k','LineWidth',2.0);
    grid on; box on;
    xlabel('热海水流量 mww [kg/s]'); ylabel('寄生功率 [kW]');
    legend({'Wp (工质泵)','Www_pump (温海水泵)','Wcw_pump (冷海水泵)','总寄生'},'Location','best');
    title('寄生功率消耗 vs 热海水流量');

    print(fig4, 'fig4_parasitic.png', '-dpng', '-r150');
    fprintf('图4已保存: fig4_parasitic.png\n');

    % ============= 综合性能 (双坐标 → 双子图) =============
    % Octave 无 yyaxis，改为并排子图
    fig5 = figure('Position',[300 300 1000 420]);
    subplot(1,2,1);
    plot(mww, R.Wnet_sys(v)/1e3, '-o','Color',C1,'LineWidth',2.0,'MarkerSize',7);
    grid on; box on;
    xlabel('热海水流量 mww [kg/s]'); ylabel('系统净发电功率 Wsys [kW]');
    title('系统净功率 vs mww');

    subplot(1,2,2);
    plot(mww, R.DT(v), '-s','Color',C2,'LineWidth',2.0,'MarkerSize',7);
    grid on; box on;
    xlabel('热海水流量 mww [kg/s]'); ylabel('总产水量 DT [kg/s]');
    title('淡水产量 vs mww');

    print(fig5, 'fig5_combined.png', '-dpng', '-r150');
    fprintf('图5已保存: fig5_combined.png\n');
end

%% =========================================================================
%       辅助函数 - 安全 LMTD
% =========================================================================
function L = safe_lmtd(dT1, dT2)
    if dT1 <= 0 || dT2 <= 0
        L = -1e3;
        return
    end
    if abs(dT1 - dT2) < 1e-8
        L = dT1;
    else
        L = (dT1 - dT2) / log(dT1/dT2);
    end
end

%% =========================================================================
%       氨 (NH3) 工质经验物性函数
% =========================================================================
function p = NH3_psat(T)
    TK = T + 273.15;
    p = 10^(4.86886 - 1113.928/(TK - 10.409)) * 100;
end

function rho = NH3_rhol(T)
    rho = 638.6 - 1.353*T - 0.00310*T^2;
end

function h = NH3_hf(T)
    h = 200.0 + 4.476*T + 0.00520*T^2;
end

function h = NH3_hg(T)
    h = 1462.0 + 1.530*T - 0.00560*T^2;
end

function s = NH3_sf(T)
    s = 1.000 + 0.01583*T - 1.80e-5*T^2;
end

function s = NH3_sg(T)
    s = 5.616 - 0.00620*T - 5.0e-6*T^2;
end

function [h3s, x] = NH3_isentropic_exp(s_in, T_out)
    sf = NH3_sf(T_out);  sg = NH3_sg(T_out);
    hf = NH3_hf(T_out);  hg = NH3_hg(T_out);
    x  = (s_in - sf) / (sg - sf);
    x  = max(0, min(1, x));
    h3s = hf + x*(hg - hf);
end
