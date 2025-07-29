%% Example 6 - Mn6
clear Sys0
rcm = 29979.2458;   meV = rcm*8.065; 

Sys0.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2

J_S4_S4   = -5*meV;    % (-5)   MnIII - MnIII.       Keep value fixed
J_S4_S5_1 = 0.4*meV;  % (0.41) MnIII - MnII, MnIII JT involved.       
J_S4_S5_2 = -0.4*meV;  % (-0.4) MnIII - MnII, MnIII JT not involved.   
J_S5_S5   = -0.1*meV; % (-0.1) MnII  - MnII.       
J_AB = J_S4_S4;
J_A1 = J_S4_S5_1; J_A3 = J_S4_S5_1; J_B2 = J_S4_S5_1; J_B4 = J_S4_S5_1; %For these, the MnIII JT axis is involved. Can be FM or AFM
J_A2 = J_S4_S5_2; J_A4 = J_S4_S5_2; J_B1 = J_S4_S5_2; J_B3 = J_S4_S5_2; %For these, the MnIII JT axis is NOT involved. Can only be AFM
J_12 = J_S5_S5;   J_34 = J_S5_S5; 
J_14 = 0;         J_23 = 0; %Assumed zero. Only interacts via a pivalate bridge.
J_13 = 0;         J_24 = 0; %Assumed zero. No interaction pathway

Sys0.J = [J_AB J_A1 J_A2 J_A3 J_A4 J_B1 J_B2 J_B3 J_B4 J_12 J_13 J_14 J_23 J_24 J_34].*(-2); %-2J formalism

DIII = -0.03*meV; %  (-0.029) MnIII anisotropy. Free Value 
EIII = 0*DIII;   % MnIII rhombicity. For INS, assume 0.
B20III = 3*DIII; B22III = EIII; %Converting to Stevens Operator formalism
DII = -0.000*meV; % MnII anisotropy. For INS, assume 0. 
EII = DII*0;      % MnII rhombicity. For INS, assume 0.
B20II = 3*DII; B22II = EII; %Converting to Stevens Operator formalism

Sys0.B2 = [B22III 0 B20III 0 0;
          B22III 0 B20III 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0];

Vary = Sys0; %Vary.J(1)=0;   %Fix value of J_AB

%Two sets of eigenvalues found - depending on the multiplicty
Exp.ev = [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.5218, 2.5218, 2.5218, 2.5218]'.*meV;   %or
Exp.ev = [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.316, 2.316,  2.5218, 2.5218]'.*meV;
Exp.ev = [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134]'.*meV;
Exp.ev = [0, 0.1414, 0.59079, 0.59079, 1.0841, 1.0841, 1.0841, 1.4134, 1.4134, 2.1868, 2.1868, 2.1868, 2.3123, 2.3123].*meV;
Exp.ev = [0, 0.1411, 0.5906, 0.5906, 1.0337, 1.0852, 1.0852, 1.4128, 1.4128, 2.1502, 2.1815, 2.1815, 2.3118, 2.3118].*meV;
%[0, 0.1411(4), 0.5906(3), 0.5906(3), 1.0337(18), 1.0852(5), 1.0852(5), 1.4128(4), 1.4128(4), 2.1502(14), 2.1815(12), 2.1815(12), 2.3118(10), 2.3118(10)]
%%
clear Opt
Opt.Verbose = true;
Opt.Scaled = true;
% Opt.C1 =1e-2;
Opt.Alpha0 = 1e0;
Opt.Linesearch = "Quadratic";
% Opt.IEPType = "Classic";
Opt.StepTolerance = 0.1;
Opt.GradientTolerance = 1e-4;
Opt.MaxIter = 200;
% Opt.Method = "RGD_LP";
Opt.NDeflations = 2 ;
% Opt.SysFound = SysFound;
SysOut = INS_IEP(Sys0,Vary,Exp,Opt)
%% Use Mint to simulate INS spectrum and compare agains baseline parameters

%Set up
Ei = 4; %Incident neutron energy in meV
MintExp.lwfwhm = 0.02*Ei/2.355; %calculated line width full-width-at-half-max, 2 % of incident energy
% MintExp.lwfwhm = 0.05;
MintExp.Energy = linspace(-Ei*0.8,Ei*0.8,1000); %calculating the spectrum in the interval -0.8*Ei to 0.8*Ei
MintExp.Q = 0.1:0.01:2.5; %Q-range which the simulation integrates over.
MintExp.Temperature = [1.5 5 10 30];

MintOpt.NumEigs = 100; %100 eigenvalues gives a good INS sim

b =[0 0.4470 0.7410];
r=[0.8500 0.3250 0.0980];
y=[0.9290 0.6940 0.1250];
g=[0.4660 0.6740 0.1880];
colours = [b;y;g;r];

for i = 1:length(SysOut)
    MintSys = SysOut(i);
    MintSys.FormFactor = {'Mn3', 'Mn3', 'Mn2', 'Mn2', 'Mn2', 'Mn2'};
    MintSys.Coords = [24.60550  13.06674   7.07523; % A
              22.27025  11.49336   6.95626; % B
              21.99544  13.57222   9.30268; % 1
              24.62949  10.95050   9.41655; % 2
              24.24486  10.55088   4.68547; % 3
              22.84040  14.03029   4.66904; % 4
              ];
    MintOpt.NumEigs = 50;
    [cross_sect] = mint(MintSys,MintExp,MintOpt);
    MintSysPrev = MintSys;
    f = figure(i);
    subplot(2,1,1)
    plot(MintExp.Energy,cross_sect*1e-4,'linewidth',1.2)
    legend('1.5 K', '5 K', '10 K', '30 K')
    xlim([0 1.5]);xlabel('E [meV]');ylabel('Signal');title('New Sim')
    ax = gca;  ax.ColorOrder = colours;
    %
    subplot(2,1,2)
    load("Sys0_Sim.mat")
    plot(MintExp.Energy,cross_sect_Sys0*1e-4,'linewidth',1.2)
    legend('1.5 K', '5 K', '10 K', '30 K')
    xlim([0 1.5]);xlabel('E [meV]');ylabel('Signal');title('Old Sim')
    ax = gca;  ax.ColorOrder = colours;
    pause
    % f.Units = 'centimeters';
    % f.Position = [10 10 20 25];
end
