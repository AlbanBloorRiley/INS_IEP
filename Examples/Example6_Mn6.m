%% Example 6 - Mn6
clear Sys0
rcm = 29979.2458;   meV = rcm*8.065; 
ExpSys.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2

% Here are some somewhat reasonable starting parameters
                 

J_S4_S4   = -5*meV;    % MnIII - MnIII.                       Rodolphe: Strong and AFM.    Keep fixed.
J_S4_S5_1 = 0.41*meV;  % MnIII - MnII, MnIII JT involved.     Rodolphe: Weak, FM or AFM.   Free value
J_S4_S5_2 = -0.4*meV; % MnIII - MnII, MnIII JT not involved. Rodolphe: Weak and AFM.      Free value
J_S5_S5   = -0.1*meV;  % MnII  - MnII.    

%Assigning parameter values to spin site pairs
%Labelling follows the drawing in the presentation from 29/4 2024
J_AB = J_S4_S4;
J_A1 = J_S4_S5_1; J_A3 = J_S4_S5_1; J_B2 = J_S4_S5_1; J_B4 = J_S4_S5_1; %For these, the MnIII JT axis is involved. Can be FM or AFM
J_A2 = J_S4_S5_2; J_A4 = J_S4_S5_2; J_B1 = J_S4_S5_2; J_B3 = J_S4_S5_2; %For these, the MnIII JT axis is NOT involved. Can only be AFM
J_12 = J_S5_S5;   J_34 = J_S5_S5; 
J_14 = 0;         J_23 = 0;% 0.01*meV; %Assumed zero. Only interacts via a pivalate bridge.
J_13 = 0;         J_24 = 0; %Assumed zero. No interaction pathway

ExpSys.J = [J_AB J_A1 J_A2 J_A3 J_A4 J_B1 J_B2 J_B3 J_B4 J_12 J_13 J_14 J_23 J_24 J_34].*(-2); %-2J formalism

DIII = -0.029*meV; % MnIII anisotropy. Free Value 
EIII = 0*DIII;   % MnIII rhombicity. For INS, assume 0.
B20III = 3*DIII; B22III = EIII; %Converting to Stevens Operator formalism

DII = -0.000*meV; % MnII anisotropy. For INS, assume 0. 
EII = DII*0;      % MnII rhombicity. For INS, assume 0.
B20II = 3*DII; B22II = EII; %Converting to Stevens Operator formalism

ExpSys.B2 = [B22III 0 B20III 0 0;
          B22III 0 B20III 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0];


FormFactor = {'Mn3', 'Mn3', 'Mn2', 'Mn2', 'Mn2', 'Mn2'};
% Metal ion coordinates for INS. Only relative positions are important.
Coords = [24.60550  13.06674   7.07523; % A
              22.27025  11.49336   6.95626; % B
              21.99544  13.57222   9.30268; % 1
              24.62949  10.95050   9.41655; % 2
              24.24486  10.55088   4.68547; % 3
              22.84040  14.03029   4.66904; % 4
              ];

% Sim INS powder spectrum
MintExp.SpectrumType = 'SE'; %INS intensity vs. energy

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



%%
Sys0.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2

J_S4_S4   = -5*meV;    % MnIII - MnIII.   

J_S4_S5_1 = 0.41*meV;  % MnIII - MnII, MnIII JT involved.     Rodolphe: Weak, FM or AFM.   Free value
J_S4_S5_2 = -0.4*meV; % MnIII - MnII, MnIII JT not involved. Rodolphe: Weak and AFM.      Free value
J_S5_S5   = -0.1*meV;  % MnII  - MnII.   


J_S4_S5_1 = 0.31*meV;  % MnIII - MnII, MnIII JT involved.     
J_S4_S5_2 = -0.3*meV; % MnIII - MnII, MnIII JT not involved. 
J_S5_S5   = -0.05*meV;  % MnII  - MnII.       
J_AB = J_S4_S4;
J_A1 = J_S4_S5_1; J_A3 = J_S4_S5_1; J_B2 = J_S4_S5_1; J_B4 = J_S4_S5_1; %For these, the MnIII JT axis is involved. Can be FM or AFM
J_A2 = J_S4_S5_2; J_A4 = J_S4_S5_2; J_B1 = J_S4_S5_2; J_B3 = J_S4_S5_2; %For these, the MnIII JT axis is NOT involved. Can only be AFM
J_12 = J_S5_S5;   J_34 = J_S5_S5; 
J_14 = 0;         J_23 = 0;% 0.01*meV; %Assumed zero. Only interacts via a pivalate bridge.
J_13 = 0;         J_24 = 0; %Assumed zero. No interaction pathway

Sys0.J = [J_AB J_A1 J_A2 J_A3 J_A4 J_B1 J_B2 J_B3 J_B4 J_12 J_13 J_14 J_23 J_24 J_34].*(-2); %-2J formalism

DIII = -0.029*meV; % MnIII anisotropy. Free Value 
DIII = -0.05*meV; % MnIII anisotropy. Free Value 
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

Vary = Sys0; Vary.J(1)=0;


Exp.ev = [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.5218, 2.5218, 2.5218, 2.5218]'.*meV;   %or
% Exp.ev = [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.316, 2.316,  2.5218, 2.5218]'.*meV;
%%
% Sys0 = ExpSys;
% Vary = Sys0;
clear Opt
Opt.Verbose = true;
% Opt.Scaled = true;
% Opt.C1 =1e-2;
Opt.Alpha0 = 1e-2;
Opt.Linesearch = "Quadratic";
% Opt.IEPType = "Classic";
% Opt.StepTolerance = 0.3;
Opt.MaxIter = 200;
% Opt.Method = "RGD_LP";
SysOut = INS_IEP(Sys0,Vary,Exp,Opt)


%% Mint Setup 

%
Sys0.Output = [];

SysOut(2) = Sys0;
for i = 1:length(SysOut)
    MintSys = SysOut(i);
    MintSys.FormFactor = FormFactor;
    MintSys.Coords = Coords;
    MintOpt.NumEigs = 50;
    [cross_sect] = mint(MintSys,MintExp,MintOpt);
    MintSysPrev = MintSys;
    figure(i)
    plotmint(cross_sect, MintExp, colours)
end
%

    % [cross_sect,Eigs,Vecs,I_nm] = mint(MintSys,MintExp,MintOpt);
    % plotmint(cross_sect, MintExp, colours)



    function plotmint(cross_sect, MintExp, colours)
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
    f.Units = 'centimeters';
    f.Position = [10 10 20 25];
    end