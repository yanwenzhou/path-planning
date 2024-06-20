function xdot = Copy_of_HRP_500t(ship,xi,ui,H,u_w,v_w)
    %% 初始化
    L   = ship.L;   % 船长（m）
    %Lpp = ship.Lpp;  % 垂线间长（m）
    B   = ship.B;   % 船宽（m）
    dA  = ship.dA;  % 艏吃水
    dF  = ship.dF;  % 艉吃水
    dm  = (dA+dF)/2;       %平均吃水
    Cb  = ship.Cb;  % 船舶方形系数：反应船体水下体型的肥胖程度 Cb = V/(L*B*d) 
    Lw  = ship.Lw;  % 设计水线长 m
    CP  = ship.CP;  % 菱形系数：表示排水体积沿船长方向的分布情况 C_p = V/(Am*L)
    %pro_num = ship.pro_num; % 螺旋桨数量 
    g   = 9.8;      % 重力加速度 m/s^2
    Am  = ship.Am;  % 中剖面面积 m^2
    HR = ship.HR;   % 舵高 m
    AR = ship.AR;   % 舵面积m^2
    DP = ship.DP;   % 螺旋桨直径m
    lam= ship.lam;  % 展弦比 lam = H_R/b_R   平均舵长和舵宽的比值
    rho= 1000;       % 水的密度(kg/m^3)
    m  = rho*ship.nabla; %质量
    PDp= ship.PDp;  % 螺距比
    P  = ship.P;
    S  = ship.S;
    delta_dotmax = 35;% 最大舵角变化率
    delta_max = 35  ; % 最大舵角指令
    np_max = 1000   ; % 最大转速指令 rpm
    TE = 1 ; % 舵机时间常数
    Tm = 5.65;
 %% 周昭明回归公式计算 附加质量和附加惯性矩《多用途贪船的操纵性预报计算》
    mx0  = 0.01*(0.398+11.97*Cb*(1+3.73*dm/B) ...
            -2.89*Cb*L/B*(1+1.13*dm/B)+0.175*Cb*(L/B)^2*(1+0.541*dm/B) ...
            -1.107*(L/B)*(dm/B))*m;    % 船体纵向附加质量 有量纲
    my0  = (0.882-0.54*Cb*(1-1.6*dm/B)-0.156*L/B*(1-0.673*Cb)+...
                0.826*(L/B)*(dm/B)*(1-0.678*(dm/B))-0.638*Cb*(L/B)*(dm/B)*...
                (1-0.669*(dm/B)))*m;   % 船体横向附加质量 有量纲
    Jzz0 = m*(L*(0.01*L*(33-76.85*Cb*(1-0.784*Cb)+3.43*(1-0.63*Cb)*L/B))/100)^2; %船体附加质量惯性矩（kg*m^2） 有量纲
    Izz = m*(0.2536*L)^2;    % 船体惯性质量矩 （kg*m^2）《On the measurement of added mass and added moment of inertia for ship motions》
   % 浅水修正附加质量和附加惯性矩
    if H/dm<3
     mx      = ((H/dm-1)^1.3+3.77+1.14*B/dm-0.233*L/dm-3.43*Cb)/(H/dm-1)^1.3*mx0;
     my      = ((H/dm-1)^0.82+0.413+0.032*B/dm+0.0129^(B/dm)^2)/(H/dm-1)^0.82*my0;
     Jzz     = ((H/dm-1)^0.82+0.413+0.0192*B/dm+0.00554^(B/dm)^2)/(H/dm-1)^0.82*Jzz0;
    else
     mx      = mx0;
     my      = my0;
     Jzz     = Jzz0;
    end
%      Izz = 5*Izz;% 修正
%      m   = 0.7*m;
    % 无量纲化附加质量
    mx_c    = mx/(rho/2*L^2*dm);  
    my_c    = my/(rho/2*L^2*dm);
    u = xi(4);        % 前进速度m/s
    v = xi(5);        % 横漂速度m/s 
    r = xi(6);        % 艏摇角速度 rad/s
    psi = xi(3);      % 船艏向 rad
    delta = xi(7);    % 舵角   rad
    np   = xi(8)/60;      % 转速 n/s
    np_order = ui(1)/60;   % 转速指令 n/s
    delta_order = ui(2);% 舵角指令 rad
    U = sqrt(u^2+v^2); % 船舶合速度
    beita = atan(-v/u); % 漂角
%   转速与舵角限制
    if abs(delta_order) >= delta_max*pi/180
            delta_order = sign(delta_order)*delta_max*pi/180;
    end
    delta_dot = (delta_order - delta)/TE;
    if abs(delta_dot) >= delta_dotmax*pi/180
        delta_dot = sign(delta_dot)*delta_dotmax*pi/180;
    end
    np_dow = np_order - np;
    if abs(np_dow) >= 40/60
       np_order = np+sign(np_dow)*40/60;
    end
    if abs(np_order) >= np_max/60
           np_order = sign(np_order)*np_max/60;
    end
    np_dot = 1/Tm*(np_order-np)*60;
 % 无量纲处理
    v_r     = v-v_w;  
    u_r     = u-u_w;      %从向量上看是减去水流的速度分量
    U_r     = sqrt(u_r^2+v_r^2);
    u_r_c = u_r/U_r;
    v_r_c = (v_r/U_r);
    r_c   = r*L/U_r;
%% 作用在船体上的流体力和力矩
    Lambda = 2*dm/L;  % λ
    tau    = (dF - dA)/dm;  % 无量纲吃水差
    lv     = Lambda/((pi/2)*Lambda +1.4*Cb*B/L);
    Cm0    = 1.11*Cb-0.07;  % 吃水差的tau=0时Cm
    Cm     = Cm0*(1+0.208*tau);   % Cm为纵轴估算图谱
    km     = 1.13902e-6;    % 淡水运动粘性系数
    Rn     = U*L/km;         % 雷诺数
    Cf     = 0.075/(log10(Rn)-2.03)^2; %  ITTC-57   计算摩擦阻力系数Cf
    z      = U/sqrt(CP*Lw); 
    Cr     = (0.1016+0.0062*z+0.2913*z^2-0.9321*z^3+2.2832*z^4-3.187*z^5+1.8895*z^6)*Am/S;  % 剩余阻力系数
    Car    = 0.0004;
    Ct0    = Cf+Cr+Car;
    Ct     = 1.3*Ct0;   % 修正

    % 浅水修正
     FH=U/sqrt(g*H);               % 浅水修正  
     if FH>0.5
        A=1.0+0.9755*FH-1.5926*FH^2;
        D=-121.875*FH+585.974*FH^2-928.4468*FH^3+490.581*FH^4;
        Ct=Ct0*(A+D*dm/H);  
     end

    Xuu  = -S*Ct/(L*dm)*(1+0.143*tau);   %贵岛考虑吃水差
    Xvr  = Cm*my_c-my_c;                     %Norrbin 回归公式
    Xvv  = 0.4*B/L-0.006*L/dm;                %松本回归公式       浅水影响不大
    Xrr  = 0.4*B/L-0.006*L/dm ;                     %松本回归公式
% 井上回归公式
    Yv   = -(pi/2*Lambda+1.4*Cb*B/L)*(1+0.67*tau^2);
%     Yv   = 1.2*Yv;
    Yr   = pi/4*Lambda*(1+0.8*tau);
    Nv   = -Lambda*(1-0.27*tau/lv);
%     Nv   = 0.88*Nv; % 修正
    Nr   = (-0.54*Lambda-Lambda^2)*(1+0.3*tau);
 % Yasuo Yoshimura 回归公式
    Yvv  = 0.048265-6.293*(1-Cb)*dm/B;
    Yrr  = 0.0045-0.445*(1-Cb)*dm/B;
    Yvr  = -0.3791+1.28*(1-Cb)*dm/B;
    Nvvr = -0.42361-3.5193*Cb*B/L+135.4668*(Cb*B/L)^2-686.3107*(Cb*B/L)^3; %改
    Nrr  = -0.0805+8.6092*(Cb*B/L)^2-36.9816*(Cb*B/L)^3;                   
    Nvrr = -0.05877+0.43958*(Cb*dm/B); 
% %     浅水修正
    if H/dm<=3
        Xvr=Xvr*(1-0.9879*dm/H+21.91123*(dm/H)^2-73.8161*(dm/H)^3+71.1409*(dm/H)^4); %修正 -2.228359002724968
        lam1e=Lambda/(dm/(2*H)*Lambda+(pi*dm/(2*H)*cot(pi*dm/(2*H)))^2.3);  
        lam1e1=Lambda/(dm/(2*H)*Lambda+(pi*dm/(2*H)*cot(pi*dm/(2*H)))^0.7);
        lam1e2=Lambda/(dm/(2*H)*Lambda+(pi*dm/(2*H)*cot(pi*dm/(2*H)))^1.7);
        Yv  = -(pi/2*lam1e+1.4*Cb*B/L)*(1+0.67*tau);       %   -1.080721443016296
        Yr  = pi/4*lam1e1*(1+0.8*tau);                     %    0.133591745389892
        Nv  = -lam1e2*(1-0.27*tau/lv);                     %   -0.307741474350297
        Nr  = -(0.54*lam1e1-lam1e1^2)*(1+0.3*tau);         %   -0.062092485649656  %         Nr  = -0.1028; 
        Yvv = Yvv*(1+14*(dm/H)^3.5);                        %   -2.704485447841801
        Yrr = Yrr*(1+3*(dm/H)^2.5);                         %   -0.068854019975778
        Yvr = Yvr*(1+3*(dm/H)^2.5);                         %   -0.772628838005657
        Nvvr= Nvvr*(1+6*(dm/H)^2.5);                        %   -0.700940677475804
        Nrr = Nrr*(1+5*(dm/H)^3.5);                         %   -0.039941018165484
        Nvrr= Nvrr*(1+6*(dm/H)^2.5);                        %    0.049332113971189
    end
    X_H = (Xuu*u_r_c^2+Xvv*v_r_c^2+Xvr*v_r_c*r_c+Xrr*r_c^2)*0.5*rho*L*dm*U_r^2;
    Y_H = (Yv*v_r_c+Yr*r_c+Yvv*abs(v_r_c)*v_r_c+Yvr*abs(v_r_c)*r_c+Yrr*abs(r_c)*r_c)*0.5*rho*L*dm*U_r^2;
    N_H = (Nv*v_r_c+Nr*r_c+Nvvr*v_r_c^2*r_c+Nvrr*v_r_c*r_c^2+Nrr*abs(r_c)*r_c)*0.5*rho*L^2*dm*U_r^2;
     
%% '''作用于螺旋桨上的力和力矩'''
    c     = -4.0;
    betap = beita-(-0.5)*r*L/U_r;  
    wP0 = 0.7*CP-0.3;      %双桨
    wP2 = wP0*exp(c*betap^2);
   
    if H/dm <= 3
        wP = 1-(1-wP2)*cos(1.7*dm/H);
    else
        wP = wP2;
    end
    wR0= 0.0229;
    Jp  = (1-wP)*u_r/(np*DP);    
    KT  = P(1)*Jp^2+P(2)*Jp+P(3);
    xR  = -0.5*L;     xB  = 0.228;
    lch = xB/L*100;  
    xH  = -(0.4+0.1*Cb)*L;   xR_c=-0.5; lR=-0.9*L;   
    tR  = 0.26118+0.0539*Cb-0.1755*Cb^2;
    gammaA = (B/dm)*(1.3*(1-Cb)-3.1*lch);
    kt     = 0.00023*(gammaA*L/DP)-0.028;
    betaR  = beita-2*xR*r/U_r;
    f      = kt*betaR;
    tP0  = 0.5*CP-0.18;      %双桨
    tP     = tP0-f;
    if H/dm <= 3
         tP =1-(1-tP)/(1-0.2*(dm/H)+0.7295*(dm/H)^2); 
         xH = xH*(1+0.3318*dm/H-3.2134*(dm/H)^2+2.5916*(dm/H)^3);  
    end
  
    X_P  = (1-tP)*rho*DP^4*np^2*KT;
    Y_P  =0;
    N_P  =0;
    %% 
    aH0  = 0.6784-1.3374*Cb+1.8891*Cb^2;
    if Jp<=0.3
        aH = aH0*Jp/0.3;
    else
        aH = aH0;
    end 
    wR    = wR0*wP/wP0;
    eps   = (1-wR)/(1-wP);
    s     = 1-u_r*(1-wP)/(np*DP*PDp);
    gamma = -22.2*(Cb*B/L)^2+0.02*(Cb*B/L)+0.68;
    fR    = 6.13*lam/(2.25+lam);
    
     if H/dm<=10
        gamma = gamma*(1+0.0161*dm/H+4.4222*(dm/H)^2-4.9825*(dm/H)^3);
        aH    = aH*(1+0.361*dm/H+1.1724*(dm/H)^2);
    end
    
    alphaR   = delta-gamma*betaR; 
 
    if delta < 0
        k2   = 0.935;
    end
    if delta >= 0
        k2   = 1.065;
    end
    
    gs   = DP/HR*0.6/eps*(2-(2-0.6/eps)*s)*s/(1-s)^2;
    UR   = U_r*(1-wR)*sqrt(1+k2*gs)/2.5;
    FN   = -0.5*rho*AR*fR*UR^2*sin(alphaR);
    X_R  = (1-tR)*FN*sin(delta);
    Y_R  = (1+aH)*FN*cos(delta);
    N_R  = (xR+aH*xH)*FN*cos(delta);
   
    X   = X_H+X_P+X_R;
    Y   = Y_H+Y_P+Y_R;
    N   = N_H+N_P+N_R;

        xdot  = [         u*cos(psi)-v*sin(psi)
                         u*sin(psi)+v*cos(psi)
                                   r
                   (X+(m+my)*v_r*r+(m+mx)*r*v_w)/(m+mx)
                    (Y-(m+mx)*u_r*r-(m+my)*r*u_w)/(m+my)
                              N/(Jzz+Izz)
                               delta_dot
                                 np_dot
                ]; 
end