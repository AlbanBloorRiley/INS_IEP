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

Sys.FormFactor = {'Mn2', 'Mn2', 'Mn2', 'Mn2'};
Sys.Coords = [21.99544  13.57222   9.30268; % 1
              24.62949  10.95050   9.41655; % 2
              24.24486  10.55088   4.68547; % 3
              22.84040  14.03029   4.66904; % 4
              ];

%% Setting up and diagonalising H

% On my machine, diagonalisation of the full H takes ~0.5 s. To describe
% the physics we probe experimentally, the first 100 eigenvalues is plenty.
% With Opt.NumEigs = 100, it takes ~0.2 s to get the needed eigenvalues.
% When computing neutron scattering, the time gained by using 'eigs' rather
% than 'eig' is even more pronounced.

clear Opt
% Opt.NumEigs = 100

tic
if exist('Opt','var') && isfield(Opt,'NumEigs')
    H = sham(Sys, [0 0 0],'sparse');
    [Vecs,E]=eigs(H,Opt.NumEigs,'smallestreal');
else
    disp('hej')
    H = sham(Sys, [0 0 0]);
    [Vecs,E]=eig(H);
end
toc
[EE, idx] = sort(real(diag(E)));
% EE = diag(E);
Vecs = Vecs(:,idx);
Eigs = (EE-min(EE))./meV; %convert to meV
Eigs(1:10)
%%
clear Exp;

Exp.SpectrumType = 'SE';
Ei = 5; %in meV, 1.28 = 8 Å
Exp.lwfwhm = 0.02*Ei/2.355;
Exp.Energy = 0:0.001:Ei*0.8;
Exp.Q = 0.1:0.01:2.5;
Exp.Temperature = [1.5 5 10 30];

if exist('Opt','var') && isfield(Opt,'NumEigs')
    [cross_sect, Eigs] = mint(Sys,Exp,Opt);
else
    [cross_sect, Eigs] = mint(Sys,Exp);
end


%%
figure
plot(Exp.Energy,cross_sect)
legend('1.5 K', '5 K', '10 K', '30 K')
