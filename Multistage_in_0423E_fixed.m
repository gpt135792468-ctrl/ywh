function [Tf,Td,Te,Tc,D,DT1,DT2,S,BPE,Tloss,ERR,ite,T1,T2,DT] = ...
    Multistage_in_0423E(N, Twwo, Tcwo, s, NEFe, NEFc, mww_in, mcw_in)
%% 输入参数意义：
% N: 级数
% Twwo: OTEC蒸发器出口温海水温度 (℃) -> 作为SLTD热源
% Tcwo: OTEC冷凝器出口冷海水温度 (℃) -> 作为SLTD冷源
% s: 进水盐度 (小数形式，如 0.035)
% NEFe/NEFc: 蒸发/冷凝非平衡温差系数 (0~1)
% mww_in: OTEC排出的温海水总流量 (kg/s)
% mcw_in: OTEC排出的冷海水总流量 (kg/s)

    n = N;
    % --- 初始化数组 ---
    mf=zeros(n+1,1); md=zeros(n+1,1); S=zeros(n+1,1);
    Tf=zeros(n+1,1); Td=zeros(n+1,1); Te=zeros(n,1); Tc=zeros(n,1);
    D=zeros(n,1); D2=zeros(n,1); DT1=zeros(n,1); DT2=zeros(n,1);
    BPE=zeros(n,1); Tloss=zeros(n,1);
    cpe=zeros(n,1); hfge=zeros(n,1);

    % --- 流量与温度边界条件 ---
    Tf(1)   = Twwo;      % 热侧进口 = OTEC温水出口
    Td(n+1) = Tcwo;      % 冷侧进口 = OTEC冷水出口
    mf(1)   = mww_in;    % 热侧总流量
    md(n+1) = mcw_in;    % 冷侧总流量
    S(1)    = s;

    ERR = 0;
    TE  = ones(n,1);

    % 初始蒸发温度线性分布
    for i = 1:n
        Te(i) = Tf(1) - i*(Tf(1) - Td(n+1) + 2)/(n+1);
        ERR = ERR + abs(1 - TE(i)/Te(i));
    end

    ite = 1;

    %% ========== 修正后的迭代主循环 ==========
    % 修改说明：
    %   1. 收敛容差从 1e-10*n 放宽到 1e-6*n (原值对本问题无法收敛)
    %   2. 最大迭代次数从 7000 增至 50000
    %   3. 引入自适应松弛因子 omega: 从 0.5 开始，随迭代缓慢衰减至最小 0.05
    %   4. 正向传递循环添加物理保护（防除零、防负值）
    % ==========================================
    while ERR > 1e-6*n && ite < 50000

        % ------ 自适应松弛因子 ------
        omega = max(0.5 * 0.9995^ite, 0.05);

        % ------ 正向传递：蒸发器（含物理保护） ------
        for i = 1:n
            Tave = (Tf(i) + Te(i)) / 2;
            BPE(i) = BPEw(S(i), Tave);

            dT_stage = Tf(i) - Te(i);
            if abs(dT_stage) < 1e-6          % 防止除零
                dT_stage = 1e-6;
            end
            DT1(i) = NEFe + BPE(i) / dT_stage;
            DT1(i) = min(DT1(i), 0.99);      % DT1 必须 < 1（物理约束）

            cpe(i)    = cpw(S(i), Tave);
            hfge(i)   = hlat(0, Tave);
            Tloss(i)  = 0.15;

            D(i) = mf(i)*cpe(i)*(1 - DT1(i))*(Tf(i) - Te(i)) / hfge(i);
            D(i) = max(D(i), 0);              % 产水不能为负

            mf(i+1) = mf(i) - D(i);
            if mf(i+1) <= 0                   % 进料不能耗尽
                mf(i+1) = 1e-6;
            end
            S(i+1)  = S(i)*mf(i)/mf(i+1);
            Tf(i+1) = Te(i) + DT1(i)*(Tf(i) - Te(i));
        end

        % ------ 保存旧 Te ------
        for i = 1:n
            TE(i) = Te(i);
        end

        % ------ 逆向传递：冷凝器 ------
        for i = 1:n
            j = (n+1) - i;
            Tc(j)   = Te(j) - Tloss(j);
            DT2(j)  = NEFc;
            Td(j)   = Tc(j) + DT2(j)*(Td(j+1) - Tc(j));

            Tave2   = (Td(j+1) + Tc(j)) / 2;
            cpc_j   = cpw(0, Tave2);
            hfgc_j  = hlat(0, Tave2);

            D2(j)   = md(j+1)*cpc_j*(1 - DT2(j))*(Tc(j) - Td(j+1)) / hfgc_j;
            md(j)   = md(j+1) + D2(j);
        end

        % ------ 更新 Te（自适应松弛） ------
        ERR = 0;
        DT  = 0;
        for i = 1:n
            Te_new = Tf(i) - 0.5*(D(i) + D2(i))*hfge(i) / (cpe(i)*mf(i)*(1 - DT1(i)));
            ERR    = ERR + abs(1 - TE(i)/Te_new);
            Te(i)  = TE(i) + omega*(Te_new - TE(i));   % 松弛更新
            Tf(i+1) = Te(i) + DT1(i)*(Tf(i) - Te(i));
            DT = DT + D(i);
        end

        ite = ite + 1;
    end

    T1 = Tf(n+1);
    T2 = Td(1);
end

%% --- 以下为原有的海水/蒸汽物性局部函数 ---

function thow = thow(S, T)
    a1=9.999e2; a2=2.034e-2; a3=-6.162e-3; a4=2.261e-5; a5=-4.657e-8;
    b1=8.02e2; b2=-2.001; b3=1.677e-2; b4=-3.060e-5; b5=-1.613e-5;
    thow = (a1+a2*T+a3*T^2+a4*T^3+a5*T^4) + ...
           (b1*S+b2*S*T+b3*S*T^2+b4*S*T^3+b5*S^2*T^2);
end

function cp = cpw(S, T)
    [T68, ~] = TEC(T);
    t  = T68 + 273.15;
    sp = S * 1000;    % 注意：S 必须是小数形式 (如 0.035)，sp = 35 g/kg
    A  =  5.328 - 9.76e-2*sp + 4.04e-4*sp^2;
    B  = -6.913e-3 + 7.351e-4*sp - 3.15e-6*sp^2;
    C  =  9.6e-6  - 1.927e-6*sp + 8.23e-9*sp^2;
    D  =  2.5e-9  + 1.66e-9*sp  - 7.125e-12*sp^2;
    cp = 1000*(A + B*t + C*t^2 + D*t^3);
end

function BPE = BPEw(S, T)
    % S 必须是小数形式 (如 0.035)
    A   = -4.584e-4*T^2 + 2.823e-1*T + 17.95;
    B   =  1.536e-4*T^2 + 5.267e-2*T +  6.56;
    BPE = A*S^2 + B*S;
end

function hfgw = hlat(S, T)
    hfg   = 2.501e6 - 2.369e3*T + 2.678e-1*T^2 - 8.013e-3*T^3 - 2.079e-5*T^4;
    hfgw  = hfg * (1 - S/1000);
end

function [T68, T48] = TEC(T)
    a  = [0 -0.148759 -0.267408 1.080760 1.269056 -4.089591 -1.871251 7.438081 -3.536296];
    DT = 0;
    for i = 1:9
        DT = DT + a(i)*(T/630)^(i-1);
    end
    T68 = T - DT;
    T48 = -113586.36365 + (113586.363652 + 227272.7273*T68)^0.5;
end
