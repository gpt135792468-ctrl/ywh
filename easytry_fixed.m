clear; clc; close all;

%% ===================================================================
%  闭式朗肯循环 (NH3) + SLTD 多级低温蒸馏 耦合系统
%  热源：OTEC温海水 (Twwi=28°C)  冷源：OTEC冷海水 (Tcwi=5°C)
%  目标：扫描热侧流量 mww，计算净功率、产水量、GOR 等性能指标
% ===================================================================

%% ===== 1. 边界条件与设备参数 =====

% --- 海水参数 ---
Twwi = 28;          % 热海水进口温度 [°C]
Tcwi =  5;          % 冷海水进口温度 [°C]
mcw  = 10;          % 冷海水流量 [kg/s]
cpw  = 4020;        % 海水比热容 [J/(kg·K)]  （参考值，全程常数）

% --- OTEC 设备参数 ---
UAe    = 50e3;      % 蒸发器 UA [W/K]
UAc    = 50e3;      % 冷凝器 UA [W/K]
eta_t  = 0.80;      % 透平等熵效率
eta_p  = 0.75;      % 工质泵效率
eta_wp = 0.80;      % 海水泵效率

% --- 工质流量 ---
mr = 0.1;           % 氨工质流量 [kg/s]

% --- 海水泵扬程 ---
% ★ 修正：Hcw 是液压扬程（压降/ρg），不是取水深度
%    800 m 深水管道的实际液压扬程约 20~30 m（管摩擦损失）
%    此处取 25 m 作为合理值；如需用 800 m，Wnet_sys 将始终为负
Hww = 10;           % 温海水泵扬程 [m]（表层取水，含管损）
Hcw = 25;           % ★ 已修正：冷海水泵液压扬程 [m]（原值 800 为取水深度，非扬程）

% --- SLTD 参数 ---
N_stage = 4;        % 级数
S0      = 0.035;    % 进水盐度（小数形式，0.035 = 35 g/kg）★ 不能写 35
NEFe    = 0.02;     % 蒸发侧非平衡温差系数
NEFc    = 0.02;     % 冷凝侧非平衡温差系数

%% ===== 2. 扫描热侧流量 =====

mww_vec = (10:5:90)';
nM      = numel(mww_vec);

% --- 预分配结果结构体 ---
fld = {'mww','Teva','Tcon','Twwo','Tcwo', ...
       'h1','h2','h3','h4','s1','s2','s3','s4', ...
       'Wt','Wp','Wnet','Wnet_sys','Qeva','Qcon','eta_cycle', ...
       'Www_pump','Wcw_pump', ...
       'DT','RR','GOR','T_sltd_hot','T_sltd_cold','converged'};
for i = 1:numel(fld)
    R.(fld{i}) = zeros(nM,1);
end

opts = optimoptions('fsolve','Display','off', ...
       'TolFun',1e-10,'TolX',1e-10, ...
       'MaxIterations',500,'MaxFunctionEvaluations',3000);

x0 = [22; 10; 25; 7];   % 初值：[Teva, Tcon, Twwo, Tcwo]

fprintf('==================================================\n');
fprintf('  闭式 OTEC + SLTD 耦合系统 - 热侧流量扫描\n');
fprintf('==================================================\n');
fprintf('开始扫描 %d 个工况 (mww = %g ~ %g kg/s)...\n\n', ...
    nM, mww_vec(1), mww_vec(end));

for k = 1:nM
    mww = mww_vec(k);

    % --- 求解 4×4 非线性方程组 ---
    [xk, fval, flag] = fsolve(@(x) residuals_LMTD(x, ...
        Twwi, Tcwi, mww, mcw, mr, UAe, UAc, cpw, eta_t, eta_p), x0, opts);

    if flag <= 0 || norm(fval) > 1e-3
        warning('mww = %.0f: fsolve 未收敛, ||F|| = %.2e', mww, norm(fval));
        R.converged(k) = 0;
        continue
    end

    Teva = xk(1);  Tcon = xk(2);
    Twwo = xk(3);  Tcwo = xk(4);

    % --- 物理合理性检查 ---
    if Teva >= Twwi || Tcon <= Tcwi || Twwo <= Teva || Tcwo >= Tcon
        warning('mww = %.0f: 温度不满足物理约束', mww);
        R.converged(k) = 0;
        continue
    end

    %% --- Rankine 循环性能 ---
    [h1,h2,h3,h4,s1,s2,s3,s4,Wt,Wp] = ...
        performance(Teva, Tcon, mr, eta_t, eta_p);

    Qeva      = mr*(h2 - h1)*1e3;          % 蒸发器热量 [W]
    Qcon      = mr*(h3 - h4)*1e3;          % 冷凝器热量 [W]
    Wnet      = Wt - Wp;                   % Rankine 净功率 [W]
    eta_cycle = Wnet / Qeva;               % 循环热效率

    g = 9.81;
    Www_pump  = mww*g*Hww / eta_wp;        % 温海水泵耗功 [W]
    Wcw_pump  = mcw*g*Hcw / eta_wp;        % 冷海水泵耗功 [W]
    Wnet_sys  = Wnet - Www_pump - Wcw_pump;% 系统净功率（扣除寄生功）[W]

    %% --- SLTD 多级蒸馏 ---
    [~,~,~,~,~,~,~,~,~,~,~,~,T1s,T2s,DT] = ...
        Multistage_in_0423E(N_stage, Twwo, Tcwo, S0, NEFe, NEFc, mww, mcw);

    % ★ 修正：Q_total 应为 SLTD 可用热量，即温水从 Twwo 降至 Tcwo 释放的热量
    %         原代码错用了 (Twwi - Twwo)（OTEC蒸发器吸热，ΔT≈3°C）
    %         正确应使用 (Twwo - Tcwo)（SLTD 驱动温差，ΔT≈17°C）
    Q_total   = mww * cpw * (Twwo - Tcwo);         % ★ 修正后

    hfg_local = 2501e3 - 2.369e3*Twwo;             % 随温度变化的汽化潜热 [J/kg]
    GOR       = DT * hfg_local / Q_total;           % 造水比 (Gain Output Ratio)
    RR        = DT / mww;                           % 回收率 (Recovery Rate) = 产水/进料

    %% --- 存储结果 ---
    R.mww(k)=mww;  R.Teva(k)=Teva;  R.Tcon(k)=Tcon;
    R.Twwo(k)=Twwo; R.Tcwo(k)=Tcwo;
    R.h1(k)=h1; R.h2(k)=h2; R.h3(k)=h3; R.h4(k)=h4;
    R.s1(k)=s1; R.s2(k)=s2; R.s3(k)=s3; R.s4(k)=s4;
    R.Wt(k)=Wt; R.Wp(k)=Wp; R.Wnet(k)=Wnet; R.Wnet_sys(k)=Wnet_sys;
    R.Qeva(k)=Qeva; R.Qcon(k)=Qcon; R.eta_cycle(k)=eta_cycle;
    R.Www_pump(k)=Www_pump; R.Wcw_pump(k)=Wcw_pump;
    R.DT(k)=DT; R.RR(k)=RR; R.GOR(k)=GOR;
    R.T_sltd_hot(k)=T1s; R.T_sltd_cold(k)=T2s;
    R.converged(k) = 1;

    x0 = xk;   % 用当前解做下一工况初值（延续法）

    if mod(k, 5) == 0
        fprintf('  进度: %d/%d (%.0f%%)\n', k, nM, 100*k/nM);
    end
end
fprintf('扫描完成!\n\n');

%% ===== 3. 打印结果表格 =====
print_table(R);

%% ===== 4. 绘制性能评估图 =====
plot_performance(R);

%% ===== 5. 保存数据 =====
% save('OTEC_SLTD_results.mat', 'R');

%% =========================================================================
%   核心函数 1 — 残差方程组（严格 LMTD）
% =========================================================================
function F = residuals_LMTD(x, Twwi, Tcwi, mww, mcw, mr, UAe, UAc, cpw, eta_t, eta_p)
    Teva = x(1);  Tcon = x(2);
    Twwo = x(3);  Tcwo = x(4);

    h2 = NH3_hg(Teva);
    s2 = NH3_sg(Teva);
    h4 = NH3_hf(Tcon);

    pe = NH3_psat(Teva);
    pc = NH3_psat(Tcon);
    vl = 1 / NH3_rhol(Tcon);
    wp = vl*(pe - pc)*1e3 / eta_p;
    h1 = h4 + wp/1000;

    [h3s, ~] = NH3_isentropic_exp(s2, Tcon);
    h3 = h2 - eta_t*(h2 - h3s);

    Qeva_r = mr*(h2 - h1)*1e3;
    Qeva_w = mww*cpw*(Twwi - Twwo);
    Qcon_r = mr*(h3 - h4)*1e3;
    Qcon_w = mcw*cpw*(Tcwo - Tcwi);

    LMTDe = safe_lmtd(Twwi - Teva, Twwo - Teva);
    LMTDc = safe_lmtd(Tcon - Tcwi, Tcon - Tcwo);

    F = zeros(4,1);
    F(1) = Qeva_r - Qeva_w;
    F(2) = Qeva_r - UAe*LMTDe;
    F(3) = Qcon_r - Qcon_w;
    F(4) = Qcon_r - UAc*LMTDc;
end

%% =========================================================================
%   核心函数 2 — 透平 / 工质泵 性能计算
% =========================================================================
function [h1,h2,h3,h4,s1,s2,s3,s4,Wt,Wp] = ...
         performance(Teva, Tcon, mr, eta_t, eta_p)

    h2 = NH3_hg(Teva);   s2 = NH3_sg(Teva);
    h4 = NH3_hf(Tcon);   s4 = NH3_sf(Tcon);

    [h3s, ~] = NH3_isentropic_exp(s2, Tcon);
    h3 = h2 - eta_t*(h2 - h3s);

    hf3 = NH3_hf(Tcon);  hg3 = NH3_hg(Tcon);
    sf3 = NH3_sf(Tcon);  sg3 = NH3_sg(Tcon);
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
%   打印表格
% =========================================================================
function print_table(R)
    fprintf('============================================================\n');
    fprintf('%5s %6s %6s %6s %6s %8s %8s %6s %7s %7s %7s\n', ...
        'mww','Teva','Tcon','Twwo','Tcwo','Wnet','Wsys','eta%','DT','RR','GOR');
    fprintf('%5s %6s %6s %6s %6s %8s %8s %6s %7s %7s %7s\n', ...
        'kg/s','°C','°C','°C','°C','kW','kW','-','kg/s','-','-');
    fprintf('------------------------------------------------------------\n');
    v = R.converged == 1;
    for k = find(v)'
        fprintf('%5.0f %6.2f %6.2f %6.2f %6.2f %8.2f %8.2f %6.3f %7.4f %7.4f %7.4f\n', ...
            R.mww(k),R.Teva(k),R.Tcon(k),R.Twwo(k),R.Tcwo(k), ...
            R.Wnet(k)/1e3, R.Wnet_sys(k)/1e3, R.eta_cycle(k)*100, ...
            R.DT(k), R.RR(k), R.GOR(k));
    end
    fprintf('============================================================\n\n');
    fprintf('说明：RR = DT/mww（回收率）；GOR = DT*hfg/Q_sltd（造水比）\n\n');
end

%% =========================================================================
%   绘图函数
% =========================================================================
function plot_performance(R)
    v   = R.converged == 1;
    mww = R.mww(v);

    C1 = [0.00 0.45 0.74];
    C2 = [0.85 0.33 0.10];
    C3 = [0.47 0.67 0.19];
    C4 = [0.49 0.18 0.56];
    C5 = [0.93 0.69 0.13];

    % ===== 主性能图 =====
    figure('Name','OTEC+SLTD 主性能图','Position',[60 60 1280 800],'Color','w');

    subplot(2,3,1);
    plot(mww, R.Wt(v)/1e3, '-o','Color',C1,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',6); hold on;
    plot(mww, R.Wnet(v)/1e3, '-s','Color',C2,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',6);
    plot(mww, R.Wnet_sys(v)/1e3, '-d','Color',C3,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',6);
    grid on; box on;
    xlabel('m_{ww} [kg/s]','FontSize',11);
    ylabel('功率 [kW]','FontSize',11);
    legend({'W_t (透平)','W_{net} (Rankine)','W_{sys} (净=透平−泵)'}, ...
        'Location','best','FontSize',9);
    title('(a) 发电功率 vs m_{ww}','FontSize',12);

    subplot(2,3,2);
    plot(mww, R.DT(v), '-o','Color',C2,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',6);
    grid on; box on;
    xlabel('m_{ww} [kg/s]','FontSize',11);
    ylabel('总产水量 D_T [kg/s]','FontSize',11);
    title('(b) SLTD 淡水产量 vs m_{ww}','FontSize',12);

    subplot(2,3,3);
    plot(mww, R.GOR(v), '-^','Color',C4,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',6);
    grid on; box on;
    xlabel('m_{ww} [kg/s]','FontSize',11);
    ylabel('GOR [-]','FontSize',11);
    title('(c) 造水比 GOR vs m_{ww}','FontSize',12);

    subplot(2,3,4);
    plot(mww, R.RR(v)*100, '-v','Color',C5,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',6);
    grid on; box on;
    xlabel('m_{ww} [kg/s]','FontSize',11);
    ylabel('回收率 RR [%]','FontSize',11);
    title('(d) 产水回收率 RR vs m_{ww}','FontSize',12);

    subplot(2,3,5);
    plot(mww, R.eta_cycle(v)*100, '-o','Color',C1,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',6);
    grid on; box on;
    xlabel('m_{ww} [kg/s]','FontSize',11);
    ylabel('循环热效率 η [%]','FontSize',11);
    title('(e) Rankine 循环效率 vs m_{ww}','FontSize',12);

    subplot(2,3,6);
    plot(mww, R.Teva(v), '-o','Color',C2,'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',5); hold on;
    plot(mww, R.Tcon(v), '-s','Color',C1,'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',5);
    plot(mww, R.Twwo(v), '--d','Color',C2,'LineWidth',1.3);
    plot(mww, R.Tcwo(v), '--^','Color',C1,'LineWidth',1.3);
    plot(mww, R.T_sltd_hot(v), ':','Color',C4,'LineWidth',1.6);
    grid on; box on;
    xlabel('m_{ww} [kg/s]','FontSize',11);
    ylabel('温度 [°C]','FontSize',11);
    legend({'T_{eva}','T_{con}','T_{wwo}','T_{cwo}','T_{SLTD,out}'}, ...
        'Location','best','FontSize',8,'NumColumns',2);
    title('(f) 关键温度 vs m_{ww}','FontSize',12);

    sgtitle('闭式 OTEC + SLTD 耦合系统 — 热侧流量敏感性分析', ...
        'FontSize',14,'FontWeight','bold');

    % ===== GOR & 产水量大图 =====
    figure('Name','GOR & 产水量','Position',[120 120 1100 450],'Color','w');

    subplot(1,2,1);
    plot(mww, R.GOR(v), '-^','Color',C4,'LineWidth',2.0,'MarkerFaceColor','w','MarkerSize',7);
    grid on; box on;
    xlabel('m_{ww} [kg/s]','FontSize',12);
    ylabel('GOR [-]','FontSize',12);
    title('造水比 GOR vs 热海水流量','FontSize',13);

    subplot(1,2,2);
    plot(mww, R.DT(v), '-o','Color',C2,'LineWidth',2.0,'MarkerFaceColor','w','MarkerSize',7);
    grid on; box on;
    xlabel('m_{ww} [kg/s]','FontSize',12);
    ylabel('总产水量 D_T [kg/s]','FontSize',12);
    title('淡水产量 vs 热海水流量','FontSize',13);

    % ===== 发电功率大图 =====
    figure('Name','发电功率','Position',[180 180 800 500],'Color','w');
    plot(mww, R.Wt(v)/1e3, '-o','Color',C1,'LineWidth',2.0,'MarkerFaceColor','w','MarkerSize',7); hold on;
    plot(mww, R.Wnet(v)/1e3, '-s','Color',C2,'LineWidth',2.0,'MarkerFaceColor','w','MarkerSize',7);
    plot(mww, R.Wnet_sys(v)/1e3, '-d','Color',C3,'LineWidth',2.0,'MarkerFaceColor','w','MarkerSize',7);
    grid on; box on;
    xlabel('m_{ww} [kg/s]','FontSize',12);
    ylabel('功率 [kW]','FontSize',12);
    legend({'W_t (透平输出)', 'W_{net} (Rankine净功率)', 'W_{sys} (系统净功率)'}, ...
        'Location','best','FontSize',10);
    title('不同热海水流量下的发电功率','FontSize',13);

    % ===== 寄生功率分解 =====
    figure('Name','寄生功率','Position',[240 240 800 500],'Color','w');
    plot(mww, R.Wp(v)/1e3, '-o','Color',C4,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',6); hold on;
    plot(mww, R.Www_pump(v)/1e3, '-s','Color',C2,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',6);
    plot(mww, R.Wcw_pump(v)/1e3, '-d','Color',C1,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',6);
    plot(mww, (R.Wp(v)+R.Www_pump(v)+R.Wcw_pump(v))/1e3, '--k','LineWidth',2.0);
    grid on; box on;
    xlabel('m_{ww} [kg/s]','FontSize',12);
    ylabel('寄生功率 [kW]','FontSize',12);
    legend({'W_p (工质泵)','W_{ww} (温海水泵)','W_{cw} (冷海水泵)','总寄生功率'}, ...
        'Location','best','FontSize',10);
    title('寄生功率消耗 vs 热海水流量','FontSize',13);

    % ===== 综合性能双坐标图 =====
    figure('Name','综合性能','Position',[300 300 800 500],'Color','w');
    yyaxis left;
    plot(mww, R.Wnet_sys(v)/1e3, '-o','Color',C1,'LineWidth',2.0,'MarkerFaceColor','w','MarkerSize',7);
    ylabel('系统净功率 W_{sys} [kW]','FontSize',12,'Color',C1);
    set(gca,'YColor',C1);
    yyaxis right;
    plot(mww, R.DT(v), '-s','Color',C2,'LineWidth',2.0,'MarkerFaceColor','w','MarkerSize',7);
    ylabel('总产水量 D_T [kg/s]','FontSize',12,'Color',C2);
    set(gca,'YColor',C2);
    grid on; box on;
    xlabel('m_{ww} [kg/s]','FontSize',12);
    title('发电功率与淡水产量综合表现','FontSize',13);
    legend({'W_{sys}','D_T'},'Location','best','FontSize',10);
end

%% =========================================================================
%   辅助函数 — 安全 LMTD
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
%   NH3 工质经验物性函数（NIST REFPROP 拟合，IIR 基准，适用 -10~40°C）
% =========================================================================

function p = NH3_psat(T)
    TK = T + 273.15;
    p  = 10^(4.86886 - 1113.928/(TK - 10.409)) * 100;
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

function hfg = NH3_hfg(T)
    hfg = NH3_hg(T) - NH3_hf(T);
end

function s = NH3_sf(T)
    s = 1.000 + 0.01583*T - 1.80e-5*T^2;
end

function s = NH3_sg(T)
    s = 5.616 - 0.00620*T - 5.0e-6*T^2;
end

function [h3s, x] = NH3_isentropic_exp(s_in, T_out)
    sf  = NH3_sf(T_out);  sg  = NH3_sg(T_out);
    hf  = NH3_hf(T_out);  hg  = NH3_hg(T_out);
    x   = (s_in - sf) / (sg - sf);
    x   = max(0, min(1, x));
    h3s = hf + x*(hg - hf);
end
