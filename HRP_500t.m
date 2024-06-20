function [X_H,Y_H,N_H,X_P,X_R,Y_R,N_R,delta_dot,n_dot] = HRP(delta_max,Delta_max,delta_c,delta,n_c,n_max,TE,ship,v_r_c,u_r_c,r_c,v_r,r,u_r,n,U,mx,my,m,beta,H,id)
    %% 初始化
    L   = ship.L;
    B   = ship.B;
    dF  = ship.dF;
    dA  = ship.dA;
    Cb  = ship.Cb;
    Lw  = ship.Lw;
    CP  = ship.CP;
    Lpp = ship.Lpp;
    pro_num = ship.pro_num;
    g   = 9.8;
    Am  = ship.Am;
    mx_c    = mx/m;
    my_c    = my/m; 
    HR = ship.HR;
    AR = ship.AR;
    DP = ship.DP;
    lam= ship.lam;
    rho= 1025;
    PDp= ship.PDp;
    P  = ship.P;
    S  = ship.S;
%     id = ship.id;
 
   %% 
    %船体小漂角采用井上模型，大漂角采用芳村模型
    %换船后需要更改的参数 S Ct
    
    d   = (dA+dF)/2;       %平均吃水
    tau = (dA-dF)/d;       %无量纲吃水差  
    lam1= 2*d/L;           %不是展弦比
    Cm0 = 1.11*Cb-0.07;    
    Cm  = Cm0*(1+0.208*tau);
    lv  = lam1/(pi/2*lam1+1.4*Cb*B/L);

    k2  = 1.13902e-6;        %流体运动粘性系数，取淡水
    Rn  = U*L/k2;            %雷诺数
     z=U/sqrt(CP*Lw);
%      Cr=(0.9006+0.6062*(z)+0.2913*(z)^2-0.2921*(z)^3+2.2832*(z)^4-2.347*(z)^5+2.5895*(z)^6)*Am/S;
     Cr=(0.1016+0.0062*(z)+0.2913*(z)^2-0.9321*(z)^3+2.2832*(z)^4-3.187*(z)^5+1.8895*(z)^6)*Am/S;
     
     Car=0.0004;
     Cf = 0.075/(log10(Rn)-2.03)^2;
     Ct0 = 1.3* (Cf+Cr+Car);       % 修正
   	 Ct =Ct0;   

     FH=U/sqrt(g*H);               % 浅水修正  
     if FH>0.5
        A=1.0+0.9755*FH-1.5926*FH^2;
        D=-121.875*FH+585.974*FH^2-928.4468*FH^3+490.581*FH^4;
        Ct=Ct0*(A+D*d/H);  
     end
    Xuu = -S*Ct/(L*d);
    Xvr = Cm*my_c-my_c;
    Xvv = 0.4*B/L-0.006*L/d; %浅水影响不大
    Xrr = -0.0003*L/d;
        
    %注意这里的Ct值        
    Yv   = -(pi/2*lam1+1.4*Cb*B/L)*(1+0.67*tau);
    Yr   = pi/4*lam1*(1+0.8*tau);
    Nv   = -lam1*(1-0.27*tau/lv);
    Nr   = -(0.54*lam1-lam1^2)*(1+0.3*tau);
    Yvv  = 0.048265-6.293*(1-Cb)*d/B;
    Yrr  = 0.0045-0.445*(1-Cb)*d/B;
    Yvr  = -0.3791+1.28*(1-Cb)*d/B;
    Nvvr = -0.42361-3.5193*Cb*B/L+135.4668*(Cb*B/L)^2-686.3107*(Cb*B/L)^3; %改
    Nrr  = -0.0805+8.6092*(Cb*B/L)^2-36.9816*(Cb*B/L)^3;                   
    Nvrr = -0.05877+0.43958*(Cb*d/B);  
     
    %% 浅水修正
    if H/d<=10
        Xvr=Xvr*(1-0.9879*d/H+21.91123*(d/H)^2-73.8161*(d/H)^3+71.1409*(d/H)^4); %修正 -2.228359002724968
        lam1e=lam1/(d/(2*H)*lam1+(pi*d/(2*H)*cot(pi*d/(2*H)))^2.3);  
        lam1e1=lam1/(d/(2*H)*lam1+(pi*d/(2*H)*cot(pi*d/(2*H)))^0.7);
        lam1e2=lam1/(d/(2*H)*lam1+(pi*d/(2*H)*cot(pi*d/(2*H)))^1.7);
        Yv  = -(pi/2*lam1e+1.4*Cb*B/L)*(1+0.67*tau);       %   -1.080721443016296
        Yr  = pi/4*lam1e1*(1+0.8*tau);                     %    0.133591745389892
        Nv  = -lam1e2*(1-0.27*tau/lv);                     %   -0.307741474350297
        Nr  = -(0.54*lam1e1-lam1e1^2)*(1+0.3*tau);         %   -0.062092485649656  %         Nr  = -0.1028; 
        Yvv = Yvv*(1+14*(d/H)^3.5);                        %   -2.704485447841801
        Yrr = Yrr*(1+3*(d/H)^2.5);                         %   -0.068854019975778
        Yvr = Yvr*(1+3*(d/H)^2.5);                         %   -0.772628838005657
        Nvvr= Nvvr*(1+6*(d/H)^2.5);                        %   -0.700940677475804
        Nrr = Nrr*(1+5*(d/H)^3.5);                         %   -0.039941018165484
        Nvrr= Nvrr*(1+6*(d/H)^2.5);                        %    0.049332113971189
    end
    
        X_H = (Xuu*u_r_c^2+Xvv*v_r_c^2+Xvr*v_r_c*r_c+Xrr*r_c^2)*0.5*rho*L*d*U^2;
        Y_H = (Yv*v_r_c+Yr*r_c+Yvv*abs(v_r_c)*v_r_c+Yvr*abs(v_r_c)*r_c+Yrr*abs(r_c)*r_c)*0.5*rho*L*d*U^2;
        N_H = (Nv*v_r_c+Nr*r_c+Nvvr*v_r_c^2*r_c+Nvrr*v_r_c*r_c^2+Nrr*abs(r_c)*r_c)*0.5*rho*L^2*d*U^2;
     
    %% Propeller  Rudder
    if abs(delta_c) >= delta_max*pi/180
        delta_c = sign(delta_c)*delta_max*pi/180;
    end
        delta_dot = (delta_c - delta)/TE;

        if abs(delta_dot) >= Delta_max*pi/180
        delta_dot = sign(delta_dot)*Delta_max*pi/180;
        end
        n_c = n_c /60;
           n_dow = n_c - n;
           if abs(n_dow) >= 40/60
              n_c = n+sign(n_dow)*40/60;
           end
        if abs(n_c) >= n_max/60
           n_c = sign(n_c)*n_max/60;
        end
%          if n > 0.3,Tm=5.65/n;else,Tm=18.83;end  
            Tm=5.65;
            
        n_dot = 1/Tm*(n_c-n)*60;
        
    c     = -4.0;
    betap = beta-(-0.5)*r*L/U;
    
    if pro_num ==1
         wP0 = 0.77*CP-0.18;    %单桨
    end
    if pro_num == 2
         wP0 = 0.7*CP-0.3;      %双桨
    end
    
    wP2 = wP0*exp(c*betap^2);
   
    if H/d <= 10
        wP = 1-(1-wP2)*cos(1.7*d/H);
    else
        wP = wP2;
    end
    wR0= 0.0229;
    Jp  = (1-wP)*u_r/(n*DP);    
    KT  = P(1)*Jp^2+P(2)*Jp+P(3);
    xR  = -0.5*L;     xB  = 0.228;
    lch = xB/L*100;  
    xH  = -(0.4+0.1*Cb)*L;   xR_c=-0.5; lR=-0.9*L;   
    tR  = 0.26118+0.0539*Cb-0.1755*Cb^2;
    gammaA = (B/d)*(1.3*(1-Cb)-3.1*lch);
    kt     = 0.00023*(gammaA*L/DP)-0.028;
    betaR  = beta-2*xR*r/U;
    f      = kt*betaR;
    
    if pro_num ==1
         tP0  = 0.5*CP-0.12;    %单桨
    end
    if pro_num == 2
         tP0  = 0.5*CP-0.18;      %双桨
    end
    tP     = tP0-f;
    if H/d <= 10
         tP =1-(1-tP)/(1-0.2*(d/H)+0.7295*(d/H)^2); 
         xH = xH*(1+0.3318*d/H-3.2134*(d/H)^2+2.5916*(d/H)^3);  
    end
  
    X_P  = (1-tP)*rho*DP^4*n^2*KT;
    if id == 2 || id ==3 ||id ==4
        X_P= 1*X_P;
    end
    %% 
    aH0  = 0.6784-1.3374*Cb+1.8891*Cb^2;
    if Jp<=0.3
        aH = aH0*Jp/0.3;
    else
        aH = aH0;
    end 
    wR    = wR0*wP/wP0;
    eps   = (1-wR)/(1-wP);
    s     = 1-u_r*(1-wP)/(n*DP*PDp);
    gamma = -22.2*(Cb*B/L)^2+0.02*(Cb*B/L)+0.68;
    fR    = 6.13*lam/(2.25+lam);
    
     if H/d<=10
        gamma = gamma*(1+0.0161*d/H+4.4222*(d/H)^2-4.9825*(d/H)^3);
        aH    = aH*(1+0.361*d/H+1.1724*(d/H)^2);
    end
    
    alphaR   = delta-gamma*betaR; 
 
    if delta < 0
        k2   = 0.935;
    end
    if delta >= 0
        k2   = 1.065;
    end
    
    gs   = DP/HR*0.6/eps*(2-(2-0.6/eps)*s)*s/(1-s)^2;
    UR   = U*(1-wR)*sqrt(1+k2*gs)/2.5;
    FN   = -0.5*rho*AR*fR*UR^2*sin(alphaR);
    if id == 2 || id ==3 ||id ==4
        FN= 2*FN;
    end
    
    X_R  = (1-tR)*FN*sin(delta);
    Y_R  = (1+aH)*FN*cos(delta);
    N_R  = (xR+aH*xH)*FN*cos(delta);
    
end

