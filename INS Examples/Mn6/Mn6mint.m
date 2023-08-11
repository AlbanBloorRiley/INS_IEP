clear all
%Define unit conversions
rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
Tesla = meV*0.116;   % Tesla         to MHz
SysOut1=[];
% [A,~,~,Ops] = SysInput(Sys);
% SysOut = Sys_Output(d,Ops,Sys.S);
% Sys.J(1) = J_S4_S4;
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
Sys2 = Sys;
Sys2.J(1) = J_S4_S4.*(-2);
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
%%
NMinima = 1;    %Minima to find
Opt = struct('NMinima',NMinima,'UseInitialGuess',1,'Method','Gauss-NewtonT1',...
    'Linesearch','Basic','Tolerance',1e-1,'p',2,'Scaled',1,...
    'minalpha',1e-4,'tau',0.1,'theta',1,'MaxIter',100,'Verbose',true);
%Now find the minima
[SysOut, NIter, Flags, Iters, FinalError]= INS_IEP(Sys,Vary,Exp,Opt)


%%
clear Exp;
% close all
NumMinima = 1;
%Loop to plot mint for all found minima
for i = 1:NumMinima
Sys1=SysOut(i);
% Sys1 = Sys
% Sys1.J(1) = J_S4_S4*(-2); %<------ Add back in the fixed variable
Sys1.FormFactor = {'Mn3', 'Mn3', 'Mn2', 'Mn2', 'Mn2', 'Mn2'};

Sys1.Coords = [24.60550  13.06674   7.07523; % A
              22.27025  11.49336   6.95626; % B
              21.99544  13.57222   9.30268; % 1
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
pause
end
%% Notes





%%


% d =
%    6.4143e+07
%        8117.6
%    1.8574e+05
%   -2.2957e+06
%     2.358e+05
%    -2.214e+06
%    7.1122e+05
 

%
% d =
% 
%    6.0675e+07
%        4831.6
%    2.7102e+05
%   -1.6438e+06
%         41801
%   -2.1781e+06
%    1.0505e+06

% 
% d =[6.0545e+07
% 4842.9
% 2.7101e+05
% -1.6399e+06
% 41916
% -2.1701e+06
% 1.0468e+06]   J_AB fixed but14/23 varried



% d=[      9581.2
%    2.7918e+05
%   -1.4512e+07
%   -1.2355e+07
%    2.1695e+07
%   -5.6097e+06]  %J_AB varried error<3e-7

% %% Calculate susceptibility using curry
% Opt.Output = 'ChiCGS'; % specify output format
% Exp.Temperature = 1:5:300;
% Exp.Field = 100; %mT
% chies = curry(Sys,Exp,Opt);
% %%
% figure()
% plot(Exp.Temperature, chies.*Exp.Temperature,'o')
% xlabel('T [K]')
% ylabel('\chi T [cm^3 K/mol]')