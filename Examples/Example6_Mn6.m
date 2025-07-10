%% Example 6 - Mn6
clear Sys0
rcm = 29979.2458;   meV = rcm*8.065; 

Sys0.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2

J_S4_S4   = -5*meV;    % (-5)
J_S4_S5_1 = 0.31*meV;  % (0.41) MnIII - MnII, MnIII JT involved.       
J_S4_S5_2 = -0.3*meV;  % (-0.4) MnIII - MnII, MnIII JT not involved.   
J_S5_S5   = -0.05*meV; % (-0.1) MnII  - MnII.       
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

%Two sets of eigenvalues found - depending on the multiplicty
Exp.ev = [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.5218, 2.5218, 2.5218, 2.5218]'.*meV;   %or
Exp.ev = [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.316, 2.316,  2.5218, 2.5218]'.*meV;
%%
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
%% Use Mint to simulate INS spectrum and compare agains baseline parameters

for i = 1:length(SysOut)
    MintSys = SysOut(i);
    MintSys.FormFactor = FormFactor;
    MintSys.Coords = Coords;
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
