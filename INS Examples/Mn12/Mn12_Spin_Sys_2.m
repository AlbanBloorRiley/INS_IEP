function Mn12_Spin_Sys_2
%conversions
rcm = 29979.2458;    % reciprocal cm to MHz
meV = rcm*8.065;     % meV           to MHz

B20 = -0.0570*meV/3; %(D = 3*B02)
B40 = (-2.78*10^-6)*meV;
B44 = (-3.2*10^-6)*meV; 
B22 = (6.8*10^-4)*meV; %(=E)

S = 10; 
Sys.S = S; 

n = (2*S+1)^1;
Sys.B2 = [B22 0 B20 0 0];        % B(k=2,q) with q = +2,+1,0,-1,-2
Sys.B4 = [B44 0 0 0 B40 0 0 0 0];  % B(k=4,q) with q = +4,+3,+2,+1,0,-1,-2,-3,-4


H = sham(Sys,[0,0,0]); 
[Vecs,E]=eig(H);
EE = diag(E);

% Exp.SpectrumType = 'SE';
% Exp.Temperature = 24;
% Exp.Energy = -1.5:0.005:1.5; Exp.lwfwhm = 0.1; Exp.Q = 0.2:0.05:1;
% Sys.FormFactor = 'Mn3';
% Sys.Coords = [0 0 0];
% [cross,Eigs] = mint(Sys,Exp);
% 
% figure
% plot(Exp.Energy,cross)


Sys.EE=EE;

A = SysInput(Sys);
A0 = sparse(n,n);

assignin('base','EE',EE)
assignin('base','A',A)
assignin('base','A0',A0)
assignin('base','Sys',Sys)