function [Sys,Vary,Exp] = Mn6_Sys 

%Define unit conversions
rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
Tesla = meV*0.116;   % Tesla         to MHz
SysOut1=[];
%% The measured eigenvalues are
exp_eigs = [0.14153, 0.59095, 0.59095, 0.9428, 0.9428, 0.9428, 1.2285+0.59095,1.2285+0.59095, 1.5128+0.9428, 2.603, 2.86];
exp_errs = [0.00016, 0.00011, 0.0003, sqrt(0.0016^2+0.00011^2), sqrt(0.0016^2+0.0003^2), 0.012, 0.02];
% exp_eigs = [0.14153, 0.59095, 0.9428, 1.2285+0.59095, 1.5128+0.9428, 2.86];


%%
% clear Sys; clear Opt
Sys.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2

% NB! Verify the number of exchange paths with/without JT axis in MnIII.

% These three should be varied in scenario 1
% Here are some somewhat reasonable starting parameters
J_S4_S4 = -5*meV;  % MnIII - MnIII. According to Rodolphe, this should be strong and AFM when JT axis is not along the MnIII-MnIII bond. Approx value 40 rcm/5 meV. This should stay fixed
% J_S4_S4 < 7.5meV

J_S4_S5_1 = 0.002*meV; % MnIII - MnII. This should have magnitude between MnII-MnII and MnIII-MnIII. Approx value a few K. MnIII JT axis involved. FM or AFM
%<0.3*meV
J_S4_S5_2 = -0.002*meV; % MnIII - MnII. This should have magnitude between MnII-MnII and MnIII-MnIII. Approx value a few K. MnIII JT axis not involved. AFM
%<0.3*meV

J_S5_S5 = 0.0625*meV; % MnII  - MnII. Rodoplhe says this should be smaller than MnII-MnIII. MnII-MnII is always AFM, approx. value ~0.5 cm-1
% 0>J_S5_S5<0.1*meV

%This parameter should stay constant in our first attempts
J_AB = J_S4_S4;        


% For now, these 8 params have the same values. Each row may need seperate values
J_A1 = J_S4_S5_1; J_A3 = J_S4_S5_1; J_B2 = J_S4_S5_1; J_B4 = J_S4_S5_1; %For these, the MnIII JT axis is involved. Can be FM or AFM
J_A2 = J_S4_S5_2; J_A4 = J_S4_S5_2; J_B1 = J_S4_S5_2; J_B3 = J_S4_S5_2; %For these, the MnIII JT axis is NOT involved. Can only be AFM

% For now, these 2 params have the same values. Each row may need seperate values
J_12 = J_S5_S5; 
J_34 = J_S5_S5; 

% For now, these parameters are zero. Each row could need its own distinct value 
J_14 = 0; 
J_23 = 0;
% J_14 = 0.1*J_S5_S5; 
% J_23 = 0.1*J_S5_S5;

% Should stay zero
J_13 = 0; J_24 = 0;

Sys.J = [J_AB J_A1 J_A2 J_A3 J_A4 J_B1 J_B2 J_B3 J_B4 J_12 J_13 J_14 J_23 J_24 J_34].*(-2); %-2JS.S formalism

% Stevens Operators for MnIII ions. As long as |J_S4_S4| >> |J_S4_S5| or 
% |J_S5_S5|, the lowest eigenvalues are indepedent of these values
% For scenario 1, these stay constant. Middle/green ones
DIII = 0*rcm; EIII = 0.1*DIII; % A notation Ramsus is more familiar with. These are fixed to zero
B20III = 3*DIII; B22III = EIII; %Converting to Stevens Operator formalism

% Stevens Operators for MnII ions. These are to be fitted. Purple ones
DII = 0.03*meV; EII = DII*1e-5; % A notation Ramsus is more familiar with. These are free parameters
B20II = 3*DII; B22II = EII; %Converting to Stevens Operator formalism

Sys.B2 = [B22III 0 B20III 0 0;
          B22III 0 B20III 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0];

%set the parameters to be varied
Vary.B2= [0 0 0 0 0;
          0 0 0 0 0;
          1 0 1 0 0;
          1 0 1 0 0;
          1 0 1 0 0;
          1 0 1 0 0];
Vary.J = [0 1 1 1 1 1 1 1 1 1 0 0 0 0 1]; 

NumEigs = 8;    %Number of experimental eigenvalues to be fit
Exp.ev = (exp_eigs(1:NumEigs)-exp_eigs(1)).*meV;%Add eigenvalues to be used to Sys
