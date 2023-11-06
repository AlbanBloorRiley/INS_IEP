function [Sys,Exp] = Mn12_Spin_Sys_3(Mn_sol,sf)%conversions
rcm = 29979.2458;    % reciprocal cm to MHz
meV = rcm*8.065;     % meV           to MHz


B20 = -0.0570*meV/3; %(D = 3*B02)
B40 = (-2.78*10^-6)*meV;
if Mn_sol == 1
    B44 = (-3.2*10^-6)*meV; 
    B22 = (6.8*10^-4)*meV; %(=E)
elseif Mn_sol == 2
    B44 = (5.1*10^-6)*meV; 
    B22 = (5.3*10^-4)*meV; %(=E)
else
    error('Please input 1 or 2, for the respective solution')
end
Sys.S = 10; 
Sys.B2 = [B22 0 B20 0 0];        % B(k=2,q) with q = +2,+1,0,-1,-2
Sys.B4 = [B44 0 0 0 B40 0 0 0 0];  % B(k=4,q) with q = +4,+3,+2,+1,0,-1,-2,-3,-4



H = ham(Sys,[0,0,0]);
[~,E]=eig(H);
EE = diag(E);
Exp.ev=EE-EE(1);
Sys.B2 = round([B22 0 B20 0 0],sf,'significant');        % B(k=2,q) with q = +2,+1,0,-1,-2
Sys.B4 = round([B44 0 0 0 B40 0 0 0 0],sf,'significant');  % B(k=4,q) with q = +4,+3,+2,+1,0,-1,-2,-3,-4

%% mint
% Exp.SpectrumType = 'SE';
% Exp.Temperature = 24;
% Exp.Energy = -1.5:0.005:1.5; Exp.lwfwhm = 0.1; Exp.Q = 0.2:0.05:1;
% Sys.FormFactor = 'Mn3';
% Sys.Coords = [0 0 0];
% %[cross,Eigs] = mint(Sys,Exp);
% %plot(Exp.Energy,cross)




