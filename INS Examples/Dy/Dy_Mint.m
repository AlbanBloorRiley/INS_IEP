% Exp.SpectrumType = 'SE';Exp.Temperature = 24;
% Exp.Energy = 0:0.05:4; Exp.lwfwhm = 0.4; Exp.Q = 0.2:0.05:1;


 rcm = 29979.2458; %rcm to MHz 
 meV = rcm*8.065; %meV to MHz
% % setup of the Hamiltonian in the usual easyspin way 


clear Sys1


% % parameters = {'B22','B20','B40','B60','Jex2','Jez2','Jex1','Jez1'};
% parameters = {'B22','B20','B40','B60','Jex2','Jex1'};
% % parameters = {'B20','B40','B60','Jex2','Jex1'};
% % parameters = {'B40','B60','Jex2','Jex1'};
% for i = 1:length(parameters)
%     X.(parameters{i}) = Y(i+1);
% end
% 
% 
% Sys1.S = [1/2 15/2 1/2];
% Sys1.B2 = [0 0 0 0 0 ; X.B22 0 X.B20 0 0 ; 0 0 0 0 0 ];
% Sys1.B4 = [0 0 0 0 0 0 0 0 0; 0 0 0 0 X.B40 0 0 0 0; 0 0 0 0 0 0 0 0 0];
% Sys1.B6 = [0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 X.B60 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0];
% Sys1.J = [X.Jex1 0 X.Jex2];
% Sys1.ee = [diag([X.Jex1,X.Jex1,X.Jez1]);zeros(3);diag([X.Jex2,X.Jex2,X.Jez2])];

% % 
Sys1 = SysOut(k);

% SE calculation
% Exp.Energy = [0:0.01:1.5]; % spectrum energy range in meV
% Exp.lwfwhm = .1; % gaussin broadening FWHM in meV
% Exp.Q = [0.1,0.5,1.0,1.5,2.0]; % Q range in A^-1
% Exp.SpectrumType = 'SE' ; % S(E) only
% Exp.Temperature = [2,9]; % Kelvin
Sys1.FormFactor = {'Cu2','Dy3','Cu2'}; %elements included i.e. Dy(3+) and neutral oxygen radical



Sys1.Coords = Dy_Coords;

Opt.NumEigs = 10;
Exp.SpectrumType = 'SE'; 
Exp.Temperature = [1, 3.5, 9]; 
Exp.Energy = -1.5:0.001:1.5; 
Exp.lwfwhm = 0.05;
Exp.Q = [0.1,0.5,1,1.5,2]; 
 [cross,Eigs] = mint(Sys1,Exp);
%f = figure;
% delete(get(f,'children'))
plot(Exp.Energy,cross);
xlim([0;Exp.Energy(end)])
xlabel('Energy (meV)') 
xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5])
ylabel('Intensity (arb. units)') 
legend('1 K','3.5 K','9 K')
