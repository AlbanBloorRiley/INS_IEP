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
Exp.ev = [0;0;0.72;0.72;1.12;1.12;1.28;1.28;1.28;1.28;1.58;1.58].*meV;
Exp.ev = [0;0;0.72;0.72;1.12;1.12;1.28;1.28;1.58;1.58].*meV;

Sys.S = [1/2 15/2 1/2];
B20 = -0.173585036568371*meV;B22 = 0.443886669455654*meV;
B40 = -0.000767621238475*meV;B60 = -0.000125378265669*meV;
Jex1 = 2*0.076588615657158*meV;Jex2 = 2*0.048947841762338*meV;

Sys.B2 = [0 0 0 0 0 ; B22 0 B20 0 0 ; 0 0 0 0 0 ];
% Sys.B2 = [0 0 0 0 0 ; 0 0 B20 0 0 ; 0 0 0 0 0 ];
Sys.B4 = [zeros(1,9); 0 0 0 0 B40 0 0 0 0; zeros(1,9)];
Sys.B6 = [zeros(1,13); 0 0 0 0 0 0 B60 0 0 0 0 0 0; zeros(1,13)];
% Sys.J = [Jex1 0 Jex2];

% Sys.ee = [diag([Jex1,Jex1,Jex1]);zeros(3);diag([Jex2,Jex2,Jex2])];
% 
Sys.ee = [diag([Jex1,Jex1,Jex1+1]);zeros(3);diag([Jex2,Jex2,Jex2+1])];
% Vary.ee = [[1,0,1;0,1,0;1,0,1];zeros(3);[1,0,1;0,1,0;1,0,1]];

%Fix the value of B22
Vary = Sys;
Vary.B2([2]) = 0;

%%
% This model of the system is rank deficien and does not converge within 
% 1000 iterations with default settings:
SysOut = INS_IEP(Sys,Vary,Exp);
%% 
%Note the warning about rank deficiency, now if we include a regularisation
%term the method converges very quickly:
clear Opt
Opt.MaxIter = 200;
Opt.Regularisation = 1e-8;
SysOut= INS_IEP(Sys,Vary,Exp,Opt);
%%
%Warnings are still given about how the matrix is singular or badly
%scaled. In this case either changing the linear solver or scaling the
%variables will solve this problem:
Opt.LinearSolver = "lsqminnorm";
SysOut= INS_IEP(Sys,Vary,Exp,Opt);
Opt = rmfield(Opt,'LinearSolver');
Opt.Scaled = true;
SysOut= INS_IEP(Sys,Vary,Exp,Opt);

%% 
%If multiple solutions are required:
clear Opt
Opt.Verbose = true;
Opt.NDeflations = 3;
Opt.Sigma = 1e-3;
% Opt.Theta = 'exp';
% Opt.SysFound = SysFound;
Opt.Scaled = true;
Opt.MaxIter = 400;
Opt.Regularisation = 1e-1;
Opt.LinearSolver = "lsqminnorm";
% Opt.Method = "Newton";
% Opt.Method = "Bad_GN";
Opt.C1 = 1e-10;
% Opt.StepTolerance = 1e-10;
% Opt.FunctionTolerance = 1e-4;
Opt.IEPType = "Classic";
[SysOut,Opt,params] = INS_IEP(Sys,Vary,Exp,Opt);
% Antiisotropic? off diagonals



%% Plotting the INS Spectrum of the minimising Spin Systems (Requires mint)
% mint can be downloaded from https://mlbakerlab.co.uk/mint/

Dy_coords =  [      0            0            0;
     -0.11166       1.9967      -1.2563;
      -1.3752      -1.3608       1.2909;
       1.4795       1.3061        1.261;
       1.5057      -1.4735       1.1273;
        1.971     0.047554      -1.2991;
     0.054185       -1.987      -1.2418;
      -1.9351     0.047795      -1.2996;
       -1.489       1.3153        1.224];

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
    plot(MintExp.Energy,cross);
    xlim([0;MintExp.Energy(end)])
    xlabel('Energy (meV)')
    xticks(0:0.1:1.5)
    ylabel('Intensity (arb. units)')
    legend('1 K','3.5 K','9 K')
    ylim([0,40]);% 1K
    pause
end



