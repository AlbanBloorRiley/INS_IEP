clear all
%Define unit conversions
rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
Tesla = meV*0.116;   % Tesla         to MHz

%%
clear Sys; clear Opt
Sys.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2

% Here are some somewhat reasonable starting parameters
J_S4_S4   = -5*meV;    % MnIII - MnIII.                       Rodolphe: Strong and AFM.    Keep fixed.
J_S4_S5_1 = 0.41*meV;  % MnIII - MnII, MnIII JT involved.     Rodolphe: Weak, FM or AFM.   Free value
J_S4_S5_2 = -0.41*meV; % MnIII - MnII, MnIII JT not involved. Rodolphe: Weak and AFM.      Free value
J_S5_S5   = -0.1*meV;  % MnII  - MnII.                        Rodoplhe: Very weak and AFM. Free value


%Assigning parameter values to spin site pairs
%Labelling follows the drawing in the presentation from 29/4 2024
J_AB = J_S4_S4;
J_A1 = J_S4_S5_1; J_A3 = J_S4_S5_1; J_B2 = J_S4_S5_1; J_B4 = J_S4_S5_1; %For these, the MnIII JT axis is involved. Can be FM or AFM
J_A2 = J_S4_S5_2; J_A4 = J_S4_S5_2; J_B1 = J_S4_S5_2; J_B3 = J_S4_S5_2; %For these, the MnIII JT axis is NOT involved. Can only be AFM
J_12 = J_S5_S5;   J_34 = J_S5_S5; 
J_14 = 0;         J_23 = 0; %Assumed zero. Only interacts via a pivalate bridge.
J_13 = 0;         J_24 = 0; %Assumed zero. No interaction pathway

Sys.J = [J_AB J_A1 J_A2 J_A3 J_A4 J_B1 J_B2 J_B3 J_B4 J_12 J_13 J_14 J_23 J_24 J_34].*(-2); %-2J formalism

DIII = -0.029*meV; % MnIII anisotropy. Free Value 
EIII = 0*DIII;   % MnIII rhombicity. For INS, assume 0.
B20III = 3*DIII; B22III = EIII; %Converting to Stevens Operator formalism

DII = -0.000*meV; % MnII anisotropy. For INS, assume 0. 
EII = DII*0;      % MnII rhombicity. For INS, assume 0.
B20II = 3*DII; B22II = EII; %Converting to Stevens Operator formalism

Sys.B2 = [B22III 0 B20III 0 0;
          B22III 0 B20III 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0];


% INS-specific information. Magnetic form factors.
Sys.FormFactor = {'Mn3', 'Mn3', 'Mn2', 'Mn2', 'Mn2', 'Mn2'};

% Metal ion coordinates for INS. Only relative positions are important.
Sys.Coords = [24.60550  13.06674   7.07523; % A
              22.27025  11.49336   6.95626; % B
              21.99544  13.57222   9.30268; % 1
              24.62949  10.95050   9.41655; % 2
              24.24486  10.55088   4.68547; % 3
              22.84040  14.03029   4.66904; % 4
              ];
%%       
Opt.NumEigs = 33; %Number of lowest eigenvalues to compute. N-fold degenerate states are counted N times.

tic %to time Hamiltonian diagonalisation
H = ham(Sys, [0 0 0],'sparse'); %setting up Hamiltonian
[Vecs,E]=eigs(H,Opt.NumEigs,'smallestreal'); %Diagonalising Hamiltonian
toc %to time Hamiltonian diagonalisation
EE = diag(E);
%Vecs = Vecs(:,idx);
Eigs = (EE-min(EE))./meV; %convert to meV and set lowest eigenvalue to 0.
Eigs % print lowest eigenvalues to command window

%% Sim INS powder spectrum
clear Exp;
clear Opt
Opt.NumEigs = 100; %100 eigenvalues gives a good INS sim

Exp.SpectrumType = 'SE'; %INS intensity vs. energy

Ei = 4; %Incident neutron energy in meV
Exp.lwfwhm = 0.02*Ei/2.355; %calculated line width full-width-at-half-max, 2 % of incident energy
Exp.Energy = linspace(-Ei*0.8,Ei*0.8,1000); %calculating the spectrum in the interval -0.8*Ei to 0.8*Ei
Exp.Q = 0.1:0.01:2.5; %Q-range which the simulation integrates over.
Exp.Temperature = [1.5 5 10 30];

[cross_sect, Eigs] = mint_ras(Sys,Exp,Opt);

%% plotting INS sim
figure
plot(Exp.Energy,cross_sect*1e-4,'linewidth',1.2)
legend('1.5 K', '5 K', '10 K', '30 K')
xlim([0 2.4])
xlabel('E [meV]')
ylabel('Signal')
title('Preliminary Sim')

%% Calculate susceptibility using curry

%NB! This code in theory calculates the susceptibility. However, the
%Hamiltonian is too large and a normal computer cannot run this snippet
clear Exp
clear Opt
Opt.Units = 'CGS';
Opt.Output = 'chimol';
Exp.Temperature = 1:1:300;
Exp.Field = 100; %mT
chi = curry(Sys,Exp,Opt);

figure()
plot(Exp.Temperature, chi.*Exp.Temperature,'o')
xlabel('T [K]')
ylabel('\chi T [cm^3 K/mol]')

figure()
plot(Exp.Temperature, chi,'o')
xlabel('T [K]')
ylabel('\chi [cm^3/mol]')