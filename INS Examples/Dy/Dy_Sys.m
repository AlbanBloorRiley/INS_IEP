function [Sys,Vary,Exp] =  Dy_Sys
rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
%Tesla = meV*0.116;   % Tesla         to MHz
clear Sys
clear A 
clear dint
Dy_Coords;
Sys.S = [1/2 15/2 1/2];

B20 = -0.173585036568371*meV;%-0.1736;%-0.197;
B22 = 0.443886669455654*meV;%0.4439;%0.066;
B40 = -0.000767621238475*meV;%3.137e-4;%1.5e-4;
B60 = -0.000125378265669*meV;%-1.8351e-6;%-1.25e-5;


Jex1 = 0.076588615657158*meV;%0.9; %Unit in meV %-0.54;
Jex2 = 0.048947841762338*meV;%0.6; %Unit in meV %-0.298;
% Sys.dint = [B20;B22;B40;B60;Jex1;Jex2];

Sys.B2 = [0 0 0 0 0 ; B22 0 B20 0 0 ; 0 0 0 0 0 ];
% Sys.B2 = [0 0 0 0 0 ; 0 0 B20 0 0 ; 0 0 0 0 0 ];
Sys.B4 = [0 0 0 0 0 0 0 0 0; 0 0 0 0 B40 0 0 0 0; 0 0 0 0 0 0 0 0 0];
Sys.B6 = [0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 B60 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0];
Sys.J = [Jex1 0 Jex2];
Vary.J = [1 0 1];
% Sys.ee = [diag([Jex1,Jex1,Jex1]);zeros(3);diag([Jex2,Jex2,Jex2])];
% 
% Sys.ee = [diag([Jex1,Jex1,Jex1+1]);zeros(3);diag([Jex2,Jex2,Jex2+1])];
% Vary.ee = [[1,0,1;0,1,0;1,0,1];zeros(3);[1,0,1;0,1,0;1,0,1]];

Vary.B2 = [zeros(1,5) ; 1 0 1 0 0 ; zeros(1,5)];
Vary.B4 = [zeros(1,9); 0 0 0 0 1 0 0 0 0; zeros(1,9)];
Vary.B6 = [zeros(1,13); 0 0 0 0 0 0 1 0 0 0 0 0 0; zeros(1,13)];

%dint = round(dint,2);
EE = [0;0;0.72;0.72;1.12;1.12;1.28;1.28].*meV;
Exp.ev=EE;
