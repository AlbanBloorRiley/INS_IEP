function [Sys, Vary, Exp] =  Dy_Sys1
rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
%Tesla = meV*0.116;   % Tesla         to MHz
clear Sys
clear A 
clear dint
Dy_Coords;
Sys.S = [15/2 1/2 1/2];
Bkq = load('Bkq.mat');
Bkq = struct2array(Bkq);
k=table2array(Bkq(:,1));
q=table2array(Bkq(:,2));
e=table2array(Bkq(:,3))';



% conversions:
rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
%Tesla = meV*0.116;   % Tesla         to MHz


%Setup of J=15/2 for Dy3+ and S=1/2 for each of the two HNN radical ligands
Sys1.S =[15/2;1/2;1/2];

Sys1.B2 = [flip(e(1:5)); zeros(1,5);zeros(1,5)].*rcm; % B(k=2,q) with q = +2,+1,0,-1,-2
Sys1.B4 = [flip(e(6:14)); zeros(1,9);zeros(1,9)].*rcm; % B(k=4,q) with q = +4,+3,+2,+1,0,-1,-2,-3,-4
Sys1.B6 = [flip(e(15:27));zeros(1,13);zeros(1,13)].*rcm; %[0 0 B64 0 0 0 B60 0 0 0 0 0 0 ; zeros(1,13); zeros(1,13)].*meV;  % B(k=6,q) with q = +6,+4,+3,+2,+1,0,-1,-2,-3,-4,-6
%Sys1.B8 = [flip(e(28:44));zeros(1,17);zeros(1,17)].*rcm;
%Sys1.B10 = [flip(e(45:65));zeros(1,21);zeros(1,21)].*rcm;
%Sys1.B12 = [flip(e(66:90));zeros(1,25);zeros(1,25)].*rcm;
N_electrons = 3;
%A0 = zfield(Sys,1:3,'sparse')
Sys1.J = zeros(1,(N_electrons-1)*(N_electrons)/2);
H = sham(Sys1,[0,0,0]); %magnetic field (mT) [x,y,z]
% [~,E]=eig(H);
% EE = diag(E);
% Eigs = (EE-EE(1))./meV; %convert to meV with the lowest eig set to zero.
A0 = H;



Jex1 = 0.076588615657158*meV;%0.9; %Unit in meV %-0.54;
Jex2 = 0.048947841762338*meV;%0.6; %Unit in meV %-0.298;
Jex1 = 2*kelvin;%0.9; %Unit in meV %-0.54;
Jex2 =  2*kelvin+10;%0.6; %Unit in meV %-0.298;

%Sys.dint = [B20;B22;B40;B60;Jex1;Jex2];

%B20=1;B40=1;B60=1;Jex1=1;Jex2=1;



% Sys.B2 = [0 0 0 0 0 ; B22 0 B20 0 0 ; 0 0 0 0 0 ];
% Sys.B2 = [0 0 0 0 0 ; 0 0 B20 0 0 ; 0 0 0 0 0 ];
% Sys.B4 = [0 0 0 0 0 0 0 0 0; 0 0 0 0 B40 0 0 0 0; 0 0 0 0 0 0 0 0 0];
% Sys.B6 = [0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 B60 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0];
%Sys.J = [Jex1 Jex2 0];
%Sys.ee = [diag([Jex1,Jex1,Jex1]);zeros(3);diag([Jex2,Jex2,Jex2])];

Sys.ee = [[Jex1,0,1;0,Jex1+2,0;-1,0,Jex1+1];diag([Jex2,Jex2+2,Jex2+1]);zeros(3)];
Sys.ee = [[Jex1,0,1;0,Jex1+2,0;-1,0,Jex1+1];[Jex2,0,2;0,Jex2+2,0;-2,0,Jex2+1];zeros(3)];


Vary.B2 = [zeros(1,5) ; 0 0 1 0 0 ; zeros(1,5)];
Vary.B4 = [zeros(1,9); 0 0 0 0 1 0 0 0 0; zeros(1,9)];
Vary.B6 = [zeros(1,13); 0 0 0 0 1 0 0 0 0; zeros(1,13)];
Vary.ee = [[1,0,1;0,1,0;1,0,1];[1,0,2;0,1,0;1,0,1];zeros(3)];
% Sys.ee = [[0,0,1;0,0,0;-1,0,0];diag([0,0,0]);zeros(3)];


%dint = round(dint,2);
EE = [0;0;0.72;0.72;1.12;1.12;1.28;1.28].*meV;
%EE = [0;0;0.72;0.72;1.12;1.12;1.28;1.28;2.4;2.4;3;3].*meV;
% Sys.ev=EE;
Exp.ev = EE;
%A0 = sparse(64,64);
% [A,dint,EE,Ops] = SysInput(Sys);


%[Y,NIter,Flags] =IEP(A,EE-EE(1),[dint],1,'GN',1e-3,1000,2)
%Output_eigs

% Good INS simulation?
  %Y =[1.4733e+07;45138;1331.5;-161.68;-12.687;50581;76467];



 %  sol_1=[7.9984e+05;-128.36;0.48175;0.26221;47906]; %without B20
%    
% 

%LM     dint to 1 s.f.
%    3.9165e+06    2.783e+06   4.2711e+06   7.1222e+06    3.571e+06
%         12121        13760        21892   1.3963e+05        15731
%       -4819.8       -32072        17416       -58227        58515
%       -6.8052       119.44      -153.52      -117.25       93.841
%     0.0039678      -1.0269      -2.0043        1.827     -0.58086
%         36127        35775  -3.8679e+05        34257        38234
%    7.4016e+05       -70834        38898        39698       -32501
%
%LM    ones  (- 9th wasnt minimum
% %     6.6909e+05   9.3491e+05    7.433e+05   2.0684e+06
%       -521.92       6626.2       32.806       1519.3
%       -1032.5      -2204.4      -1165.3      -2264.2
%       -14.224      -27.267      0.70026      -6.5508
%      -0.01538      0.05889      0.20204     0.021997
%        -59026        40330        52378  -4.2724e+05
%        -24647       -39002       -39454        34754
% 
%   Columns 5 through 8
% 
%    8.1034e+05   1.3899e+06   6.4539e+05   7.3081e+05
%       -6043.3      -1403.1      -208.44       3772.8
%        2398.9       4820.3      -563.51       1517.8
%       -19.544      -74.045      -12.468      -3.8323
%     -0.045858      -0.3232     -0.17233      0.15077
%         23121       -41728       -23723       -36863
%         37718       -62772        35471       -53737
% 
%   Columns 9 through 10
% 
%        1.1367e+06
%        -2510.6
%        -1105.1
%        -18.907
%        0.022817
%        -43866
%        -1.5094e+05