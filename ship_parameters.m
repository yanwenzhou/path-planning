 function ship = ship_parameters(id)
    switch id
        case 1
      %       2000吨干散货船参数
              ship.L=67.6;
              ship.Lpp=62;   %垂线间长
              ship.Lw=120;   %设计水线长。待定
              ship.B=13.8;
              ship.dF=3.0;
              ship.dA=3.3;
              ship.CP=0.734; 
              ship.Cb=0.7;   %对回转半径影响显著
              ship.AR=7.06;   %13.5 舵面积
              ship.HR=3.28;   % 舵高
              ship.DP=1.57;
              ship.lam=1.68;
              ship.nabla=2000;
              ship.pro_num=2;  %双桨
              ship.PDp=0.68;   %螺距比
              ship.AeAo=0.44;  %盘面比
              ship.Am=28.152;  %中剖面面积
              ship.P = [ -0.0694   -0.4441    0.5277];
              ship.S = 1487;
              ship.id=1;
        case 2
      %       1000吨干散货船参数 设计航速14.56km/h- 4.04m/s
              ship.L=55.95;
              ship.Lpp=55.0; %垂线间长
              ship.Lw=100;    %设计水线长。待定
              ship.B=9.8;
              ship.dF=2.85;
              ship.dA=2.85;
              ship.CP=0.734; 
              ship.Cb=0.806;   %对回转半径影响显著
              ship.AR=4.212;    %13.5
              ship.HR=2.48;
              ship.DP=1.082;
              ship.lam=1.68;
              ship.nabla=1140;
              ship.pro_num=2;  %双桨
              ship.PDp=0.90;   %螺距比
              ship.AeAo=0.55;  %盘面比
              ship.Am=28.152;  %中剖面面积
              ship.P = [ -0.0694   -0.4441    0.5277];
              ship.S = 859.44;
              ship.id=2;
        case 3
              % 48TEU集装箱船  1000r/min 3.75m/s
              ship.L=88;%54.7
              ship.Lpp=88; %垂线间长
              ship.Lw=88;    %设计水线长。待定
              ship.B=10.8; % 12.1
              ship.dF=2; %2.55
              ship.dA=2; % 2.55
              ship.CP=0.839; 
              ship.Cb=0.834;   %对回转半径影响显著
              ship.AR=3.5;    %13.5
              ship.HR=2;
              ship.DP=1.508;
              ship.lam=1.68;
              ship.nabla=1401.2;
              ship.pro_num=2;  %双桨
              ship.PDp=0.683;   %螺距比
              ship.AeAo=0.454;  %盘面比
              ship.Am=30.67;  %中剖面面积
              ship.P = [ -0.0694   -0.1141    0.2891];
              ship.S = 1392.8;
              ship.id=3;
 %       1000吨集装箱船 4.2m/s
%               ship.L=49.34;
%               ship.Lpp=48.14; %垂线间长
%               ship.Lw=110;    %设计水线长。待定
%               ship.B=12.5;
%               ship.dF=2.48;
%               ship.dA=2.48;
%               ship.CP=0.784; 
%               ship.Cb=0.6898;   %对回转半径影响显著
%               ship.AR=4.212;    %13.5
%               ship.HR=2.48;
%               ship.DP=1.736;
%               ship.lam=1.68;
%               ship.nabla=1078;
%               ship.pro_num=2;  %双桨
%               ship.PDp=0.673;   %螺距比
%               ship.AeAo=0.55;  %盘面比
%               ship.Am=32.152;  %中剖面面积
%               ship.P = [ -0.0694   -0.4441    0.5277];
%               ship.S = 812.44;
%               ship.id=3;
        case 4
%       500吨散货船 16.3km/h 4.5 拖轮
              ship.L=111;
              ship.Lpp=44.6; %垂线间长
              ship.Lw=111;    %设计水线长。待定
              ship.B=10.8;
              ship.dF=1.6;
              ship.dA=1.6;
              ship.CP=0.694; 
              ship.Cb=0.755;   %对回转半径影响显著
              ship.AR=3.212;    %13.5
              ship.HR=1.68;
              ship.DP=0.9444;
              ship.lam=1.68;
              ship.nabla=1300;
              ship.pro_num=2;  %双桨
              ship.PDp=0.673;   %螺距比
              ship.AeAo=0.465;  %盘面比
              ship.Am=14.762;  %中剖面面积
               ship.P = [-0.219,-0.204,0.275];
%               ship.P = [ -0.0694   -0.4441    0.5277];
              ship.S = 812.44;
              ship.id=4;
        case 5
%       唐河500吨散货船 16.3km/h 4.5
              ship.L=57.5;
              ship.Lpp=51; %垂线间长
              ship.Lw=100;    %设计水线长。待定
              ship.B=10.8;
              ship.dF=1.6;
              ship.dA=1.6;
              ship.CP=0.694; 
              ship.Cb=0.649;   %对回转半径影响显著
              ship.AR=3.212;    %13.5
              ship.HR=1.68;
              ship.DP=0.9444;
              ship.lam=1.68;
              ship.nabla=650;
              ship.pro_num=2;  %双桨
              ship.PDp=0.673;   %螺距比
              ship.AeAo=0.465;  %盘面比
              ship.Am=11.9232;  %中剖面面积
               ship.P = [-0.219,-0.204,0.575];
%               ship.P = [ -0.0694   -0.4441    0.5277];
              ship.S = 845.44;
              ship.id=5;

         case 6
      %       唐河 2x500t 驳船 111.0×10.8×1.6
              ship.L=111;
              ship.Lpp=111; %垂线间长
              ship.Lw=111;    %设计水线长。待定
              ship.B=10.8;
              ship.dF=1.6;
              ship.dA=1.6;
              ship.CP=0.694; 
              ship.Cb=0.649;   %对回转半径影响显著
              ship.AR=3.212;    %13.5
              ship.HR=1.68;
              ship.DP=0.9444;
              ship.lam=1.68;
              ship.nabla=1000;
              ship.pro_num=2;  %双桨
              ship.PDp=0.673;   %螺距比
              ship.AeAo=0.465;  %盘面比
              ship.Am=11.9232;  %中剖面面积
               ship.P = [-0.219,-0.204,0.575];
%               ship.P = [ -0.0694   -0.4441    0.5277];
              ship.S = 845.44;
              ship.id=6;
    end

end

