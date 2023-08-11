% Exp.SpectrumType = 'SE';Exp.Temperature = 24;
% Exp.Energy = 0:0.05:4; Exp.lwfwhm = 0.4; Exp.Q = 0.2:0.05:1;


 rcm = 29979.2458; %rcm to MHz meV = rcm*8.065; %meV to MHz
% % setup of the Hamiltonian in the usual easyspin way 
% D = -0.041.*meV; % convert from meV to MHz
% E = 0.007.*meV; % convert from meV to MHz
E=Y(3);
D=3*Y(2);
J=Y(4);
 Sys1.S = [3/2 3/2 3/2 3/2 3/2 3/2];
 Sys1.D = [D E;D E;D E;D E;D E;D E];
 Sys1.J = [J 0 0 0 0 J 0 0 0 J 0 0 J 0 J]; % 
% 



Sys1.FormFactor = 'Cr3';



a = [0;-3.3807;-5.7290;-5.5960;-3.1148;0.2482];    
b = [0;0.0840;-2.3170;-5.6470;-7.8270;-7.7410];
c = [0; -0.2765 ; 0.1103 ; -0.2956 ; 0.2920 ; -0.1803];

Dy_Coords
Sys1.Coords = Dy1;

Opt.NumEigs = 24;
Exp.SpectrumType = 'SE'; 
Exp.Temperature = [1.5, 6, 15]; 
Exp.Energy = -1:0.001:4; 
Exp.lwfwhm = 0.05;
Exp.Q = [0.1,0.5,1,1.5,2]; [cross,Eigs] = mint(Sys1,Exp,Opt);
figure
plot(Exp.Energy,cross)
xlim([0 3.25])
xlabel('Energy (meV)') 
ylabel('Intensity (arb. units)') 
legend('1.5 K','6 K','15 K')


% 
% 
% [cross,Eigs] = mint(Sys1,Exp);
% 
% figure
% plot(Exp.Energy,cross)