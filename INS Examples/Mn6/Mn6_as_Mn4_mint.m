clear all
%Define unit conversions
rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
Tesla = meV*0.116;   % Tesla         to MHz

%% The measured eigenvalues are
exp_eigs = [0.14153, 0.59095, 0.9428, 1.2285+0.59095, 1.5128+0.9428, 2.603, 2.86];
exp_errs = [0.00016, 0.00011, 0.0003, sqrt(0.0016^2+0.00011^2), sqrt(0.0016^2+0.0003^2), 0.012, 0.02];


% The eigenvalues are determined by fitting gaussians to model the expe-
% rimental data using a chi^2 method. When the eigenvalues are measured 
% directly, the µ of the fitted gaussian is taken eigenvalues, and the
% standard error of the fit parameter is taken as the uncertainty of the
% eigenvalue. If an eigenvalue was measured in more than one experimental
% setting, the value obtained from fitting to the measurement with the best 
% resolution was taken.
% When the eigenvalue was measued indirectly, the sum of the two fitted µ's
% leading observed state is taken as the eigenvalue. In that case, the
% eigenvalue error is estimated by adding the standard error on the fitting
% parameters in quadrature. Again, the 

%% Defining system with reasonable starting params
clear Sys;
Sys.S = [5/2 5/2 5/2 5/2]; %MnII has S = 5/2

% These two J's should be varied in scenario 1
J_same = -0.3*meV;  % Interactions between MnIIs on the same side of the diamagnetic core
J_diff = -0.1001*meV; % Interactions between MnIIs on opposite sides of the diamagnetic core 

J12 = J_same; J13 = 0; J14 = J_diff; J23 = J_diff; J24 = 0; J34 = J_same;
Sys.J = [J12 J13 J14 J23 J24 J34].*(-2); %-2JS.S formalism

% Stevens Operators for MnII ions. These are to be fitted
DII = -0.05*rcm; EII = 0.1*DII; % A notation Ramsus is more familiar with
B20II = 3*DII; B22II = EII; %Converting to Stevens Operator formalism

Sys.B2 = [B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0];



NumEigs = 7;    %Number of experimental eigenvalues to be fit
Sys.ev = exp_eigs(1:NumEigs)-exp_eigs(1);%Add eigenvalues to be used to Sys
Sys.ev = (exp_eigs(1:NumEigs)-exp_eigs(1)).*meV;%Add eigenvalues to be used to Sys

NMinima = 3;    %Minima to find

%Now find the minima
[SysOut, NIter, Flags, Iters, FinalError]= INS_IEP2(Sys, 'NMinima',NMinima,'UseInitialGuess',true,'StepMethod','GN','Linesearch','No','SysFound',[],'Deflate','F','Tolerance',1e-1,'p',2);


% %% Setting up and diagonalising H
% if exist('Opt','var') && isfield(Opt,'NumEigs')
%     H = sham(Sys, [0 0 0],'sparse');
%     [Vecs,E]=eigs(H,Opt.NumEigs,'smallestreal');
% else
%     disp('hej')
%     H = sham(Sys, [0 0 0]);
%     [Vecs,E]=eig(H);
% end

%%
clear Exp;
% close all
NumMinima = NMinima;
%Loop to plot mint for all found minima
for i = 1:NumMinima
Sys1=SysOut(i);
Sys1.FormFactor = {'Mn2', 'Mn2', 'Mn2', 'Mn2'};
Sys1.Coords = [21.99544  13.57222   9.30268; % 1
              24.62949  10.95050   9.41655; % 2
              24.24486  10.55088   4.68547; % 3
              22.84040  14.03029   4.66904; % 4
              ];

Exp.SpectrumType = 'SE';
Ei = 5; %in meV, 1.28 = 8 Å
Exp.lwfwhm = 0.02*Ei/2.355;
Exp.Energy = 0:0.001:Ei*0.8;
Exp.Q = 0.1:0.01:2.5;
Exp.Temperature = [1.5 5 10 30];

 Opt.NumEigs = 50;
if exist('Opt','var') && isfield(Opt,'NumEigs')
    [cross_sect] = mint(Sys1,Exp,Opt);
else
    [cross_sect] = mint(Sys1,Exp);
end



figure
plot(Exp.Energy,cross_sect)
legend('1.5 K', '5 K', '10 K', '30 K')
% pause
end

