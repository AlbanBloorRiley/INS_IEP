clear all
%Define unit conversions
rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
Tesla = meV*0.116;   % Tesla         to MHz

%%
clear Sys; clear Opt
Sys.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2

% NB! Verify the number of exchange paths with/without JT axis in MnIII.

% These three should be varied in scenario 1
% Here are some somewhat reasonable starting parameters
J_S4_S4 = -5*meV;  % MnIII - MnIII. According to Rodolphe, this should be strong and AFM when JT axis is not along the MnIII-MnIII bond. Approx value 40 rcm/5 meV. This should stay fixed

J_S4_S5_1 = 0.2001*meV; % MnIII - MnII. This should have magnitude between MnII-MnII and MnIII-MnIII. Approx value a few K. MnIII JT axis involved. FM or AFM
J_S4_S5_2 = -0.2001*meV; % MnIII - MnII. This should have magnitude between MnII-MnII and MnIII-MnIII. Approx value a few K. MnIII JT axis not involved. AFM

J_S5_S5 = -0.0625*meV; % MnII  - MnII. Rodoplhe says this should be smaller than MnII-MnIII. MnII-MnII is always AFM, approx. value ~0.5 cm-1

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

% Should stay zero
J_13 = 0; J_24 = 0;

Sys.J = [J_AB J_A1 J_A2 J_A3 J_A4 J_B1 J_B2 J_B3 J_B4 J_12 J_13 J_14 J_23 J_24 J_34].*(-2); %-2JS.S formalism

% Stevens Operators for MnIII ions. As long as |J_S4_S4| >> |J_S4_S5| or 
% |J_S5_S5|, the lowest eigenvalues are indepedent of these values
% For scenario 1, these stay constant. Middle/green ones
DIII = 0*rcm; EIII = 0.1*DIII; % A notation Ramsus is more familiar with. These are fixed to zero
B20III = 3*DIII; B22III = EIII; %Converting to Stevens Operator formalism

% Stevens Operators for MnII ions. These are to be fitted. Purple ones
DII = 0.0125*meV; EII = DII*0; % A notation Ramsus is more familiar with. These are free parameters
B20II = 3*DII; B22II = EII; %Converting to Stevens Operator formalism

Sys.B2 = [B22III 0 B20III 0 0;
          B22III 0 B20III 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0];


A     = [25.87028  11.29944   7.52276]-[23.46093  14.95621   6.63648]; %diff between "top" and "bottom" oxygen on the JT axis for A
B     = [23.44574   9.61166   6.60192]-[20.96949  13.30406   7.29219]; %diff between "top" and "bottom" oxygen on the JT axis for B

% direct = (a+b).*units/2; direct = direct./norm(direct); %average direction of the two JT axes

% Calculation of rotation matrix rotating 'direct' into [0; 0; 1]
% v = cross(direct,[0;0;1]);
% vmat = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
% theta = asin(norm(cross(direct,[0;0;1])));
% c = dot(direct,[0;0;1])*cos(theta);
% R = eye(3)+vmat+vmat*vmat./(1+c);


% eulerangles = eulang(R); %calculates Euler angles from a rotation matrix

% Each from specifies the orientation of the local anisotropy axes in the
% molecular frame. Locally, the anisotropy axis is [0; 0; 1]
% Sys.DFrame = [eulerangles; eulerangles; 0 0 0; 0 0 0; 0 0 0; 0 0 0];


Sys.FormFactor = {'Mn3', 'Mn3', 'Mn2', 'Mn2', 'Mn2', 'Mn2'};

Sys.Coords = [24.60550  13.06674   7.07523; % A
              22.27025  11.49336   6.95626; % B
              21.99544  13.57222   9.30268; % 1
              24.62949  10.95050   9.41655; % 2
              24.24486  10.55088   4.68547; % 3
              22.84040  14.03029   4.66904; % 4
              ];
%%
% I found that the issues with 'eigs' was happened because I had not added
% the option 'smallestreal'
NumEigs = 50;
Opt.NumEigs = NumEigs; 

tic
if exist('Opt','var') && isfield(Opt,'NumEigs')
    H = sham(Sys, [0 0 0],'sparse');
    [Vecs,E]=eigs(H,Opt.NumEigs,'smallestreal');
else
    H = sham(Sys, [0 0 0]);
    [Vecs,E]=eig(H);
end
toc
[EE, idx] = sort(real(diag(E)));
% EE = diag(E);
Vecs = Vecs(:,idx);
Eigs = (EE-min(EE))./meV; %convert to meV
Eigs
%%
clear Exp;
clear Opt
NumEigs = 100;
Opt.NumEigs = NumEigs;
Exp.SpectrumType = 'SE';
Ei = 5; %in meV, 1.28 = 8 Å
Exp.lwfwhm = 0.02*Ei/2.355;
Exp.Energy = 0:0.001:Ei*0.8;
Exp.Q = 0.1:0.01:2.5;
Exp.Temperature = [1.5];

[cross_sect, Eigs] = mint_ras(Sys,Exp,Opt);

%%
figure
plot(Exp.Energy,cross_sect)
legend('1.5 K', '5 K', '10 K', '30 K')

%% Calculate susceptibility using curry
Opt.Output = 'ChiCGS'; % specify output format
Exp.Temperature = 1:5:300;
Exp.Field = 100; %mT
chies = curry(Sys,Exp,Opt);
%%
figure()
plot(Exp.Temperature, chies.*Exp.Temperature,'o')
xlabel('T [K]')
ylabel('\chi T [cm^3 K/mol]')