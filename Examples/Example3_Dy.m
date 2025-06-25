%% Dy_6 Dysprosium
clear Sys
% Simulation and Comparison  6A energy
%Only need to run this once!
% cmf = mfilename('fullpath')
% cf = fileparts(cmf)
% folder = fullfile(cf,'Dy_data')
%
%    folder = '/Volumes/mlbakerlab/neutron_data/FRMII/6396_TOFTOF/analysis/DyHNN2_analysis/2021/data/feb2022/for_Alan/';
%    folder = '/Users/user/Dropbox (The University of Manchester)/alban.bloorriley@manchester.ac.uk’s files/MatLab/mint/Dy_data/';
%    folder = '/Dy_data/';
   filename = 'DyHNN2_IofE_6A_1K_2theta_gt_30_bkgdsub=0p1.inx';
%    Dy_6A1K = importdata([folder,filename]);
   Dy_6A1K = importdata(filename);
   filename = 'DyHNN2_IofE_6A_9K_2theta_gt_30_bkgdsub_0p1.inx';
%    Dy_6A9K = importdata([folder,filename]);
    Dy_6A9K = importdata(filename);
   A= {Dy_6A1K,Dy_6A9K};
   %for i=1:length(A)
   n=2;
       x_avg_1 = arrayfun(@(j) mean(A{1}(j:j+n-1,1)),1:n:length(A{1}(:,1))-n+1)';
    y_avg_1 = arrayfun(@(j) mean(A{1}(j:j+n-1,2)),1:n:length(A{1}(:,1))-n+1)';
    err_avg_1 = arrayfun(@(j) sqrt(mean(A{1}(j:j+n-1,3).^2)),1:n:length(A{1}(:,1))-n+1)'; 
    x_avg_2 = arrayfun(@(j) mean(A{2}(j:j+n-1,1)),1:n:length(A{2}(:,1))-n+1)';
    y_avg_2 = arrayfun(@(j) mean(A{2}(j:j+n-1,2)),1:n:length(A{2}(:,1))-n+1)';
    err_avg_2 = arrayfun(@(j) sqrt(mean(A{2}(j:j+n-1,3).^2)),1:n:length(A{2}(:,1))-n+1)'; 
   %end


rcm = 29979.2458;   meV = rcm*8.065; 
%Experimental eigenvalues
Exp.ev = [0;0;0.72;0.72;1.12;1.12;1.28;1.28].*meV;

Sys.S = [1/2 15/2 1/2];
B20 = -0.173585036568371*meV;B22 = 0.443886669455654*meV;
B40 = -0.000767621238475*meV;B60 = -0.000125378265669*meV;
Jex1 = 2*0.076588615657158*meV;Jex2 = 2*0.048947841762338*meV;
% 
% B20 = -0.173585*meV;B22 =  0.443886*meV;
% B40 = -0.000767*meV;B60 = -0.000125*meV;
% Jex1 = 2*0.076588*meV;Jex2 = 2*0.048947*meV;

% B20 = -0.1*meV;B22 =  0.1*meV;
% B40 = -0.001*meV;
% B60 = -0.0001*meV;
% Jex1 = 2*0.1*meV;Jex2 = 2*0.01*meV;
Sys.B2 = [0 0 0 0 0 ; B22 0 B20 0 0 ; 0 0 0 0 0 ];
% Sys.B2 = [0 0 0 0 0 ; 0 0 B20 0 0 ; 0 0 0 0 0 ];
Sys.B4 = [zeros(1,9); 0 0 0 0 B40 0 0 0 0; zeros(1,9)];
Sys.B6 = [zeros(1,13); 0 0 0 0 0 0 B60 0 0 0 0 0 0; zeros(1,13)];
% Sys.J = [Jex1 0 Jex2];

% Sys.ee = [diag([Jex1,Jex1,Jex1]);zeros(3);diag([Jex2,Jex2,Jex2])];
% 
Sys.ee = [diag([Jex1,Jex1,Jex1+1]);zeros(3);diag([Jex2,Jex2,Jex2+1])];
% Vary.ee = [[1,0,1;0,1,0;1,0,1];zeros(3);[1,0,1;0,1,0;1,0,1]];
Vary = Sys;
%
Vary.B2([2,8]) = 0;
% Vary.B4(14) = 0;
Vary.B6(20) = 0;
%%
% This model of the system is rank deficien and does not converge within 
% 1000 iterations with default settings:
SysOut1 = INS_IEP(Sys,Vary,Exp);
%Note the warning about rank deficiency, now if we include a regularisation
%term the method converges very quickly:

%%
clear Opt
Opt.Verbose = true;
Opt.NDeflations = 15;
Opt.Sigma = 1e-5;
% Opt.Theta = 'exp';
% Opt.SysFound = SysFound;
Opt.Scaled = true;
Opt.MaxIter = 200;
% Opt.Regularisation = 1e-8;
Opt.LinearSolver = "lsqminnorm";
Opt.Method = "Newton";
% Opt.Method = "RGD_LP";
Opt.Epsilon = 1e-4;
% Opt.C1 = 1e-10;
Opt.StepTolerance = 1e-10;
Opt.FunctionTolerance = 1e-4;
% Opt.IEPType = "Classic";
[SysOut,Opt,params] = INS_IEP(Sys,Vary,Exp,Opt);
% Antiisotropic? off diagonals



%%
   Dy_coords = Dy_Coords;

   rcm = 29979.2458;    % reciprocal cm to
   
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
%Tesla = meV*0.116;   % Tesla         to MHz


MintOpt.NumEigs = 20;
MintExp.SpectrumType = 'SE'; 
MintExp.Temperature = [1, 3.5, 9]; 
MintExp.Energy = -1.5:0.001:1.5; 
MintExp.lwfwhm = 0.05;
MintExp.Q = [0.1,0.5,1,1.5,2]; 

for i = 1:length(SysOut)
    if ~contains(params.convergence.ConvergenceFlags,SysOut(i).Output.ConvergenceFlag)
        continue
    end
    SysOut(i).FormFactor = {'Cu2','Dy3','Cu2'}; %elements included i.e. Dy(3+) and neutral oxygen radical
    SysOut(i).Coords = Dy_coords;
    [cross,Eigs] = mint(SysOut(i),MintExp,MintOpt);
    i
    %f = figure;
    % delete(get(f,'children'))
    plot(MintExp.Energy,cross);
    xlim([0;MintExp.Energy(end)])
    xlabel('Energy (meV)')
    xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5])
    ylabel('Intensity (arb. units)')
    legend('1 K','3.5 K','9 K')

    hold on
    errorbar(x_avg_1, y_avg_1, err_avg_1,'o-','LineWidth',1.5);
    errorbar(x_avg_2, y_avg_2, err_avg_2,'o-','LineWidth',1.5);
    %
    xlim([0,1.5]); %1K
    ylim([0,80]);% 1K
    legend('Sim 1 K','Sim 3.5 K','Sim 9 K','Exp 1K','Exp 9K')
    hold off
    pause
end



%% New easyspin
MintExp.Field = 100;
MintExp.Temperature = Dy_Susceptibility(:,1);
[~,chizz] = curry(SysOut(1),MintExp);

%one of these
figure(1)
plot(MintExp.Temperature, chizz)
figure(2)
plot(MintExp.Temperature, MintExp.Temperature'.*chizz./(4*pi*1e-6))

%plot Dy_Susceptibility for the exp data

















%%
clear all
%Calculate eigenvalues
rcm = 29979.2458;   meV = rcm*8.065;
ExpSysA.S = 6;
ExpSysA.B2 = [0 0 0.705 0 0.0250].*meV;
ExpSysA.B4 = [0 0 0 0 0 0 0 0 -8.5E-4 ].*meV;
EE = eig(ham(ExpSysA,[0,0,0]));
ExpA.ev = EE-EE(1);
ExpSysB.S = 6;
ExpSysB.B2 = [0 0 0.712 0 0.0729].*meV;
EE = eig(ham(ExpSysB,[0,0,0]));
ExpB.ev = EE(1:end)-EE(1);

%Set up system fo 
Sys0.S = 6;
Sys0.B2 = [0 0 1 0 0.1].*meV;
% Sys0.B2 = ExpSysB.B2;
Sys0.B2 = [0 0 0.71 0 0.072].*meV;
%In this case the standard method fails to converge this is most likely due to the highly non
% linear nature of the problem:
Opt.Verbose = true;

SysOutB = INS_IEP(Sys0,Sys0,ExpB)

%In this case it is better to use an alternaitve minimisation method, 
%for example the Riemannian Gradient Descent Lift and Projection method [6],
% which finds the minimum for each example:
Opt = struct('Method', 'RGD_LP');
SysOutB = INS_IEP(Sys0,Sys0,ExpB,Opt)

Sys0.B4 = [zeros(1,8),-1e-2].*meV;
SysOutA = INS_IEP(Sys0,Sys0,ExpA,Opt)


