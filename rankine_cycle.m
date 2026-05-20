%% ============================================================
%  朗肯循环热力学模型  (Rankine Cycle Thermodynamic Model)
%  ============================================================
%
%  状态点 / State Points:
%    1 - 冷凝器出口 / 泵入口     — 饱和液体
%    2 - 泵出口 / 锅炉入口       — 压缩液体
%    3 - 锅炉出口 / 汽轮机入口   — 过热蒸汽
%    4 - 汽轮机出口 / 冷凝器入口 — 湿蒸汽或过热蒸汽
%
%  依赖 / Requires: XSteam (MATLAB 水蒸气热物性库)
%    下载: https://www.mathworks.com/matlabcentral/fileexchange/9817
%    下载后将 XSteam.m 放入 MATLAB 工作路径即可使用
%
%  压力单位约定: 对外参数用 kPa，调用 XSteam 时转换为 bar (÷100)
% ============================================================

clear; clc; close all;

%% ====== 1. 输入参数 / Input Parameters ======

P_cond_kPa = 10;       % 冷凝器压力 [kPa]
P_boil_kPa = 3000;     % 锅炉压力   [kPa]
T3_C       = 400;      % 汽轮机入口温度 [°C]  (需高于饱和温度)
eta_pump   = 0.85;     % 泵等熵效率   [-]
eta_turb   = 0.85;     % 汽轮机等熵效率 [-]
m_dot      = 1;        % 质量流量     [kg/s]

% 转换为 bar 供 XSteam 使用
Pc_bar = P_cond_kPa / 100;
Pb_bar = P_boil_kPa / 100;

%% ====== 2. 状态点计算 / State Point Calculations ======

% ---------- 状态 1: 饱和液体 (冷凝器出口 / 泵入口) ----------
T1 = XSteam('Tsat_p', Pc_bar);      % 饱和温度        [°C]
h1 = XSteam('hL_p',   Pc_bar);      % 饱和液比焓      [kJ/kg]
s1 = XSteam('sL_p',   Pc_bar);      % 饱和液比熵      [kJ/kg·K]
v1 = XSteam('vL_p',   Pc_bar);      % 饱和液比体积    [m³/kg]

% ---------- 状态 2: 压缩液体 (泵出口 / 锅炉入口) ----------
% 等熵泵功: w_p,s = v1·ΔP  (m³/kg × kPa = kJ/kg)
wp_s = v1 * (P_boil_kPa - P_cond_kPa);  % 等熵泵功  [kJ/kg]
wp   = wp_s / eta_pump;                   % 实际泵功  [kJ/kg]
h2   = h1 + wp;
T2   = XSteam('T_ph', Pb_bar, h2);
s2   = XSteam('s_ph', Pb_bar, h2);

% ---------- 状态 3: 过热蒸汽 (锅炉出口 / 汽轮机入口) ----------
T3 = T3_C;
h3 = XSteam('h_pT', Pb_bar, T3);
s3 = XSteam('s_pT', Pb_bar, T3);

% 检查是否超热
T_sat_boil = XSteam('Tsat_p', Pb_bar);
if T3 <= T_sat_boil
    error('汽轮机入口温度 T3=%.1f°C 未超过锅炉饱和温度 %.1f°C，请提高 T3_C。', T3, T_sat_boil);
end

% ---------- 状态 4: 汽轮机出口 (冷凝器入口) ----------
s4s  = s3;                                 % 等熵假设: s4s = s3
h4s  = XSteam('h_ps', Pc_bar, s4s);       % 等熵出口比焓  [kJ/kg]
wt   = eta_turb * (h3 - h4s);             % 实际汽轮机功  [kJ/kg]
h4   = h3 - wt;
T4   = XSteam('T_ph', Pc_bar, h4);
s4   = XSteam('s_ph', Pc_bar, h4);
x4   = XSteam('x_ph', Pc_bar, h4);        % 干度 (两相区); 否则 NaN or <0

%% ====== 3. 循环性能 / Cycle Performance ======

Q_in   = h3 - h2;                % 锅炉吸热     [kJ/kg]
Q_out  = h4 - h1;                % 冷凝器散热   [kJ/kg]
W_net  = wt - wp;                % 循环净功     [kJ/kg]
eta_th = W_net / Q_in * 100;     % 热效率       [%]
BWR    = wp / wt * 100;          % 回功比       [%]

% 卡诺效率 (理想上限参考)
T_H_K     = T3 + 273.15;         % 高温热源温度 [K]
T_L_K     = T1 + 273.15;         % 低温热源温度 [K]
eta_Carnot = (1 - T_L_K / T_H_K) * 100;

%% ====== 4. 输出结果 / Print Results ======

sep = repmat('=', 1, 60);
fprintf('\n%s\n', sep);
fprintf('        朗肯循环分析结果 / Rankine Cycle Results\n');
fprintf('%s\n\n', sep);

fprintf('【运行参数 / Operating Parameters】\n');
fprintf('  锅炉压力        : %8.1f kPa  (%.2f MPa)\n', P_boil_kPa, P_boil_kPa/1e3);
fprintf('  冷凝压力        : %8.1f kPa\n', P_cond_kPa);
fprintf('  汽轮机入口温度  : %8.1f °C\n', T3);
fprintf('  泵等熵效率      : %8.0f%%\n', eta_pump*100);
fprintf('  汽轮机等熵效率  : %8.0f%%\n\n', eta_turb*100);

fprintf('【状态点参数 / State Point Properties】\n');
fprintf('  %-4s  %-12s  %-14s  %-14s  %s\n', ...
        '点', 'T [°C]', 'h [kJ/kg]', 's [kJ/kg·K]', '状态');
fprintf('  %s\n', repmat('-', 1, 60));
fprintf('  %-4s  %-12.3f  %-14.3f  %-14.5f  %s\n', '1', T1, h1, s1, '饱和液体');
fprintf('  %-4s  %-12.3f  %-14.3f  %-14.5f  %s\n', '2', T2, h2, s2, '压缩液体');
fprintf('  %-4s  %-12.3f  %-14.3f  %-14.5f  %s\n', '3', T3, h3, s3, '过热蒸汽');
if isnan(x4) || x4 < 0 || x4 > 1
    state4_str = '过热蒸汽';
else
    state4_str = sprintf('湿蒸汽  x=%.4f', x4);
end
fprintf('  %-4s  %-12.3f  %-14.3f  %-14.5f  %s\n', '4', T4, h4, s4, state4_str);

fprintf('\n【能量分析 / Energy Analysis (per kg)】\n');
fprintf('  汽轮机输出功    : %8.2f kJ/kg\n', wt);
fprintf('  泵消耗功        : %8.2f kJ/kg\n', wp);
fprintf('  循环净功        : %8.2f kJ/kg\n', W_net);
fprintf('  锅炉吸热量      : %8.2f kJ/kg\n', Q_in);
fprintf('  冷凝器散热量    : %8.2f kJ/kg\n', Q_out);
fprintf('  能量守恒校验    : Q_in - Q_out - W_net = %.4f kJ/kg\n', Q_in - Q_out - W_net);

fprintf('\n【循环性能 / Cycle Performance】\n');
fprintf('  热效率 η_th     : %8.2f%%\n', eta_th);
fprintf('  卡诺效率 (参考) : %8.2f%%\n', eta_Carnot);
fprintf('  回功比 BWR      : %8.2f%%\n', BWR);

fprintf('\n【功率输出 / Power Output (ṁ = %.1f kg/s)】\n', m_dot);
fprintf('  汽轮机功率      : %8.2f kW\n', wt * m_dot);
fprintf('  泵功率          : %8.2f kW\n', wp * m_dot);
fprintf('  净功率          : %8.2f kW\n', W_net * m_dot);
fprintf('  锅炉热功率      : %8.2f kW\n', Q_in * m_dot);
fprintf('%s\n\n', sep);

%% ====== 5. 辅助数据: 饱和曲线 / Saturation Dome Data ======

P_dome = logspace(log10(0.006), log10(220.6), 400);  % bar, 三相点→临界点
T_dome  = arrayfun(@(p) XSteam('Tsat_p', p), P_dome);
sL_dome = arrayfun(@(p) XSteam('sL_p',   p), P_dome);
sV_dome = arrayfun(@(p) XSteam('sV_p',   p), P_dome);
hL_dome = arrayfun(@(p) XSteam('hL_p',   p), P_dome);
hV_dome = arrayfun(@(p) XSteam('hV_p',   p), P_dome);
ok = ~isnan(T_dome) & ~isnan(sL_dome) & ~isnan(sV_dome);

% 锅炉等压过程 (2→3): 逐步计算 T(h) 和 s(h)
h_boil = linspace(h2, h3, 120);
T_boil = arrayfun(@(h) XSteam('T_ph', Pb_bar, h), h_boil);
s_boil = arrayfun(@(h) XSteam('s_ph', Pb_bar, h), h_boil);

% 冷凝器等压过程 (4→1)
h_cond = linspace(h4, h1, 60);
T_cond = arrayfun(@(h) XSteam('T_ph', Pc_bar, h), h_cond);
s_cond = arrayfun(@(h) XSteam('s_ph', Pc_bar, h), h_cond);

%% ====== 6. T-s 图 / T-s Diagram ======

fig1 = figure('Name', '朗肯循环 T-s 图', ...
              'Position', [50 50 760 530], 'Color', 'w');

hold on; grid on; box on;

% 饱和区填充与边界
fill([sL_dome(ok), fliplr(sV_dome(ok))], ...
     [T_dome(ok),  fliplr(T_dome(ok))], ...
     [0.75 0.88 0.97], 'EdgeColor', 'none', 'FaceAlpha', 0.45);
hDome = plot([sL_dome(ok), fliplr(sV_dome(ok))], ...
             [T_dome(ok),  fliplr(T_dome(ok))], ...
             '-', 'Color', [0.1 0.3 0.75], 'LineWidth', 1.8);

% 过程线
hPump = plot([s1, s2], [T1, T2], 'k-',  'LineWidth', 2.5);           % 泵 1→2
hBoil = plot(s_boil,    T_boil,  'r-',  'LineWidth', 2.5);            % 锅炉 2→3
hTurb = plot([s3, s4],  [T3, T4],'Color',[0.0 0.60 0.15],'LineWidth',2.5); % 汽轮机 3→4
hCond = plot(s_cond,    T_cond,  'm-',  'LineWidth', 2.5);            % 冷凝器 4→1

% 等熵虚线 (理想汽轮机参考)
hIso  = plot([s3, s3], [T3, XSteam('T_ph',Pc_bar,h4s)], ...
             '--', 'Color',[0.5 0.5 0.5], 'LineWidth', 1.2);

% 状态点标记
state_s = [s1, s2, s3, s4];
state_T = [T1, T2, T3, T4];
plot(state_s, state_T, 'o', 'MarkerSize', 9, ...
     'MarkerFaceColor', [0.92 0.15 0.15], 'MarkerEdgeColor', 'k', 'LineWidth', 1.3);

% 状态点标签
off_s = [-0.15, -0.15,  0.06,  0.06];
off_T = [-18,    15,    10,   -20 ];
for i = 1:4
    text(state_s(i)+off_s(i), state_T(i)+off_T(i), num2str(i), ...
         'FontSize', 13, 'FontWeight', 'bold');
end

xlabel('比熵  s  [kJ/(kg·K)]', 'FontSize', 12);
ylabel('温度  T  [°C]',        'FontSize', 12);
title({'\bf朗肯循环 T-s 图', ...
       sprintf('\\eta_{th}=%.1f%%   P_{boil}=%.0f kPa   T_3=%.0f°C   \\eta_p=%.0f%%   \\eta_t=%.0f%%', ...
               eta_th, P_boil_kPa, T3, eta_pump*100, eta_turb*100)}, ...
      'FontSize', 12);
legend([hDome, hPump, hBoil, hTurb, hCond, hIso], ...
       {'饱和曲线','泵 1→2','锅炉 2→3','汽轮机 3→4','冷凝器 4→1','等熵参考线'}, ...
       'Location','northwest','FontSize',10);
ylim([min(state_T)-35, max(state_T)+50]);
set(gca,'FontSize',11);

%% ====== 7. h-s 图 (莫里尔图) / Mollier Diagram ======

fig2 = figure('Name', '朗肯循环 h-s 图 (莫里尔图)', ...
              'Position', [830 50 760 530], 'Color', 'w');

hold on; grid on; box on;

fill([sL_dome(ok), fliplr(sV_dome(ok))], ...
     [hL_dome(ok), fliplr(hV_dome(ok))], ...
     [0.75 0.88 0.97], 'EdgeColor', 'none', 'FaceAlpha', 0.45);
plot(sL_dome(ok), hL_dome(ok), '-', 'Color',[0.1 0.3 0.75], 'LineWidth', 1.8);
plot(sV_dome(ok), hV_dome(ok), '-', 'Color',[0.1 0.3 0.75], 'LineWidth', 1.8);

% 冷凝器过程 h-s
h_cond2 = linspace(h4, h1, 60);
s_cond2 = arrayfun(@(h) XSteam('s_ph', Pc_bar, h), h_cond2);

plot([s1, s2],  [h1, h2],  'k-', 'LineWidth', 2.5);           % 泵
plot(s_boil,    h_boil,    'r-', 'LineWidth', 2.5);            % 锅炉
plot([s3, s4],  [h3, h4],  'Color',[0.0 0.60 0.15],'LineWidth',2.5); % 汽轮机
plot(s_cond2,   h_cond2,   'm-', 'LineWidth', 2.5);            % 冷凝器
% 等熵汽轮机参考线
plot([s3, s3], [h3, h4s], '--', 'Color',[0.5 0.5 0.5], 'LineWidth', 1.2);

state_h = [h1, h2, h3, h4];
plot(state_s, state_h, 'o', 'MarkerSize', 9, ...
     'MarkerFaceColor', [0.92 0.15 0.15], 'MarkerEdgeColor', 'k', 'LineWidth', 1.3);

off_h = [-70, 50, 35, -80];
for i = 1:4
    text(state_s(i)+off_s(i), state_h(i)+off_h(i), num2str(i), ...
         'FontSize', 13, 'FontWeight', 'bold');
end

xlabel('比熵  s  [kJ/(kg·K)]', 'FontSize', 12);
ylabel('比焓  h  [kJ/kg]',     'FontSize', 12);
title('\bf朗肯循环 h-s 图 (莫里尔图)', 'FontSize', 13);
legend({'饱和区','饱和液线','饱和汽线', ...
        '泵 1→2','锅炉 2→3','汽轮机 3→4','冷凝器 4→1','等熵参考线','状态点'}, ...
       'Location','northwest','FontSize',10);
set(gca,'FontSize',11);

%% ====== 8. 参数分析 / Parametric Study ======

fig3 = figure('Name', '朗肯循环参数分析', ...
              'Position', [50 620 1100 420], 'Color', 'w');

% ---- 子图1: 热效率 vs 锅炉压力 ----
subplot(1, 2, 1);
P_arr = linspace(200, 15000, 80);
eta_arr = nan(size(P_arr));
for i = 1:numel(P_arr)
    Pb_i = P_arr(i)/100;
    try
        if XSteam('Tsat_p', Pb_i) >= T3_C; continue; end
        h1_ = XSteam('hL_p', Pc_bar);
        v1_ = XSteam('vL_p', Pc_bar);
        h2_ = h1_ + v1_*(P_arr(i) - P_cond_kPa)/eta_pump;
        h3_ = XSteam('h_pT', Pb_i, T3_C);
        s3_ = XSteam('s_pT', Pb_i, T3_C);
        h4s_= XSteam('h_ps', Pc_bar, s3_);
        if isnan(h4s_) || h4s_ <= 0; continue; end
        h4_ = h3_ - eta_turb*(h3_ - h4s_);
        eta_arr(i) = (h3_ - h4_ - (h2_ - h1_)) / (h3_ - h2_) * 100;
    catch; end
end
plot(P_arr/1e3, eta_arr, 'b-', 'LineWidth', 2); hold on;
plot(P_boil_kPa/1e3, eta_th, 'ro', 'MarkerSize', 10, ...
     'MarkerFaceColor','r', 'DisplayName','当前设计点');
xlabel('锅炉压力 [MPa]', 'FontSize', 11);
ylabel('热效率 η_{th} [%]', 'FontSize', 11);
title('热效率 vs 锅炉压力', 'FontSize', 12);
grid on; legend({'η_{th}','当前设计点'}, 'Location','southeast','FontSize',10);
set(gca,'FontSize',10);

% ---- 子图2: 热效率 vs 汽轮机入口温度 ----
subplot(1, 2, 2);
T_arr = linspace(150, 700, 80);
eta_arr2 = nan(size(T_arr));
for i = 1:numel(T_arr)
    try
        if T_arr(i) <= T_sat_boil; continue; end
        h3_ = XSteam('h_pT', Pb_bar, T_arr(i));
        s3_ = XSteam('s_pT', Pb_bar, T_arr(i));
        h4s_= XSteam('h_ps', Pc_bar, s3_);
        if isnan(h4s_) || h4s_ <= 0; continue; end
        h4_ = h3_ - eta_turb*(h3_ - h4s_);
        eta_arr2(i) = (h3_ - h4_ - wp) / (h3_ - h2) * 100;
    catch; end
end
plot(T_arr, eta_arr2, 'r-', 'LineWidth', 2); hold on;
plot(T3, eta_th, 'bo', 'MarkerSize', 10, 'MarkerFaceColor','b');
xlabel('汽轮机入口温度 T_3 [°C]', 'FontSize', 11);
ylabel('热效率 η_{th} [%]', 'FontSize', 11);
title('热效率 vs 汽轮机入口温度', 'FontSize', 12);
grid on; legend({'η_{th}','当前设计点'}, 'Location','southeast','FontSize',10);
set(gca,'FontSize',10);

sgtitle('朗肯循环参数分析 / Parametric Study', ...
        'FontSize', 14, 'FontWeight', 'bold');

fprintf('所有图表已生成完毕。\n');
