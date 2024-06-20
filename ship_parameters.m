 function ship = ship_parameters(id)
    switch id
        case 1
      %       2000�ָ�ɢ��������
              ship.L=67.6;
              ship.Lpp=62;   %���߼䳤
              ship.Lw=120;   %���ˮ�߳�������
              ship.B=13.8;
              ship.dF=3.0;
              ship.dA=3.3;
              ship.CP=0.734; 
              ship.Cb=0.7;   %�Ի�ת�뾶Ӱ������
              ship.AR=7.06;   %13.5 �����
              ship.HR=3.28;   % ���
              ship.DP=1.57;
              ship.lam=1.68;
              ship.nabla=2000;
              ship.pro_num=2;  %˫��
              ship.PDp=0.68;   %�ݾ��
              ship.AeAo=0.44;  %�����
              ship.Am=28.152;  %���������
              ship.P = [ -0.0694   -0.4441    0.5277];
              ship.S = 1487;
              ship.id=1;
        case 2
      %       1000�ָ�ɢ�������� ��ƺ���14.56km/h- 4.04m/s
              ship.L=55.95;
              ship.Lpp=55.0; %���߼䳤
              ship.Lw=100;    %���ˮ�߳�������
              ship.B=9.8;
              ship.dF=2.85;
              ship.dA=2.85;
              ship.CP=0.734; 
              ship.Cb=0.806;   %�Ի�ת�뾶Ӱ������
              ship.AR=4.212;    %13.5
              ship.HR=2.48;
              ship.DP=1.082;
              ship.lam=1.68;
              ship.nabla=1140;
              ship.pro_num=2;  %˫��
              ship.PDp=0.90;   %�ݾ��
              ship.AeAo=0.55;  %�����
              ship.Am=28.152;  %���������
              ship.P = [ -0.0694   -0.4441    0.5277];
              ship.S = 859.44;
              ship.id=2;
        case 3
              % 48TEU��װ�䴬  1000r/min 3.75m/s
              ship.L=88;%54.7
              ship.Lpp=88; %���߼䳤
              ship.Lw=88;    %���ˮ�߳�������
              ship.B=10.8; % 12.1
              ship.dF=2; %2.55
              ship.dA=2; % 2.55
              ship.CP=0.839; 
              ship.Cb=0.834;   %�Ի�ת�뾶Ӱ������
              ship.AR=3.5;    %13.5
              ship.HR=2;
              ship.DP=1.508;
              ship.lam=1.68;
              ship.nabla=1401.2;
              ship.pro_num=2;  %˫��
              ship.PDp=0.683;   %�ݾ��
              ship.AeAo=0.454;  %�����
              ship.Am=30.67;  %���������
              ship.P = [ -0.0694   -0.1141    0.2891];
              ship.S = 1392.8;
              ship.id=3;
 %       1000�ּ�װ�䴬 4.2m/s
%               ship.L=49.34;
%               ship.Lpp=48.14; %���߼䳤
%               ship.Lw=110;    %���ˮ�߳�������
%               ship.B=12.5;
%               ship.dF=2.48;
%               ship.dA=2.48;
%               ship.CP=0.784; 
%               ship.Cb=0.6898;   %�Ի�ת�뾶Ӱ������
%               ship.AR=4.212;    %13.5
%               ship.HR=2.48;
%               ship.DP=1.736;
%               ship.lam=1.68;
%               ship.nabla=1078;
%               ship.pro_num=2;  %˫��
%               ship.PDp=0.673;   %�ݾ��
%               ship.AeAo=0.55;  %�����
%               ship.Am=32.152;  %���������
%               ship.P = [ -0.0694   -0.4441    0.5277];
%               ship.S = 812.44;
%               ship.id=3;
        case 4
%       500��ɢ���� 16.3km/h 4.5 ����
              ship.L=111;
              ship.Lpp=44.6; %���߼䳤
              ship.Lw=111;    %���ˮ�߳�������
              ship.B=10.8;
              ship.dF=1.6;
              ship.dA=1.6;
              ship.CP=0.694; 
              ship.Cb=0.755;   %�Ի�ת�뾶Ӱ������
              ship.AR=3.212;    %13.5
              ship.HR=1.68;
              ship.DP=0.9444;
              ship.lam=1.68;
              ship.nabla=1300;
              ship.pro_num=2;  %˫��
              ship.PDp=0.673;   %�ݾ��
              ship.AeAo=0.465;  %�����
              ship.Am=14.762;  %���������
               ship.P = [-0.219,-0.204,0.275];
%               ship.P = [ -0.0694   -0.4441    0.5277];
              ship.S = 812.44;
              ship.id=4;
        case 5
%       �ƺ�500��ɢ���� 16.3km/h 4.5
              ship.L=57.5;
              ship.Lpp=51; %���߼䳤
              ship.Lw=100;    %���ˮ�߳�������
              ship.B=10.8;
              ship.dF=1.6;
              ship.dA=1.6;
              ship.CP=0.694; 
              ship.Cb=0.649;   %�Ի�ת�뾶Ӱ������
              ship.AR=3.212;    %13.5
              ship.HR=1.68;
              ship.DP=0.9444;
              ship.lam=1.68;
              ship.nabla=650;
              ship.pro_num=2;  %˫��
              ship.PDp=0.673;   %�ݾ��
              ship.AeAo=0.465;  %�����
              ship.Am=11.9232;  %���������
               ship.P = [-0.219,-0.204,0.575];
%               ship.P = [ -0.0694   -0.4441    0.5277];
              ship.S = 845.44;
              ship.id=5;

         case 6
      %       �ƺ� 2x500t ���� 111.0��10.8��1.6
              ship.L=111;
              ship.Lpp=111; %���߼䳤
              ship.Lw=111;    %���ˮ�߳�������
              ship.B=10.8;
              ship.dF=1.6;
              ship.dA=1.6;
              ship.CP=0.694; 
              ship.Cb=0.649;   %�Ի�ת�뾶Ӱ������
              ship.AR=3.212;    %13.5
              ship.HR=1.68;
              ship.DP=0.9444;
              ship.lam=1.68;
              ship.nabla=1000;
              ship.pro_num=2;  %˫��
              ship.PDp=0.673;   %�ݾ��
              ship.AeAo=0.465;  %�����
              ship.Am=11.9232;  %���������
               ship.P = [-0.219,-0.204,0.575];
%               ship.P = [ -0.0694   -0.4441    0.5277];
              ship.S = 845.44;
              ship.id=6;
    end

end

