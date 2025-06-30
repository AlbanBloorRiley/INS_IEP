%% Example 3 _an exchange coupled Tb – Cu single molecule magnet including anisotropic spin exchange and nuclear hyperfine coupling.
% PHYSICAL REVIEW B 88, 064405 (2013)
clear all
% some constants
rcm = 29979.2458; % reciprocal cm to MHz
meV = rcm*8.065; % meV to MHz
% All parameters come from PHYSICAL REVIEW B 88, 064405 (2013)
%Setup of J=6 for Tb3+ and S=1/2 for Cu
%% Simulate Eigenvalues:
ExpSys.S = [6 1/2];
% Exchange coupling
% Jxx = -0.54;
% Jyy = -0.54;
% Jzz = -0.298;
% ExpSys.ee = [Jxx Jyy Jzz].*meV;
ExpSys.ee = [-0.54 -0.54 -0.298].*meV;
%Stevens parameters
% B22 = 0.066;    B20 = -0.197;  
% B44 = 5e-4;     B40 = 1.5e-4;  
% B64 = 1.17e-4;  B60 = -1.25e-5;
% ExpSys.B2 = [B22 0 B20 0 0 ; zeros(1,5)].*meV;
% ExpSys.B4 = [B44 0 0 0 B40 0 0 0 0 ; zeros(1,9)].*meV;
% ExpSys.B6 = [0 0 1.17e-4 0 0 0 B60 0 0 0 0 0 0 ; zeros(1,13)].*meV; 
ExpSys.B2 = [0.066 0 -0.197 0 0 ; zeros(1,5)].*meV;
ExpSys.B4 = [5e-4 0 0 0 1.5e-4 0 0 0 0 ; zeros(1,9)].*meV;
ExpSys.B6 = [0 0 B64 0 0 0 -1.25e-5 0 0 0 0 0 0 ; zeros(1,13)].*meV; 
H = ham(ExpSys,[0,0,0]);
EE = eig(H);
%only save a subset of the eigenvalues
% Exp.ev = EE(1:10)-EE(1);
Exp.ev = EE(1:end)-EE(1);

% ExpSys.Nucs = '159Tb,63Cu';
% ExpSys.A = [6.2e-3 6.2e-3 6.2e-3 0 0 0 ; 0 0 0 0.043 0.043 0.043].*meV;
% H2 = ham(ExpSys,[0,0,0]);
%% Set up the model spin system and initial values
% Sys0.S = [6 1/2];
% Sys0.ee = [-0.5 -0.5 -0.25].*meV;
% Sys0.B2 = [0.1 0 -0.1 0 0 ; zeros(1,5)].*meV;
% Sys0.B4 = [1e-4 0 0 0 1e-4 0 0 0 0 ; zeros(1,9)].*meV;
% Sys0.B6 = [0 0 1e-4 0 0 0 -1e-4 0 0 0 0 0 0 ; zeros(1,13)].*meV; 

Sys0.S = [6 1/2];
Sys0.ee = [-0.5 -0.5 -0.05].*meV;
Sys0.B2 = [0.1 0 -0.1 0 0 ; zeros(1,5)].*meV;
Sys0.B4 = [1e-10 0 0 0 1e-10 0 0 0 0 ; zeros(1,9)].*meV;
Sys0.B6 = [0 0 1e-10 0 0 0 -1e-10 0 0 0 0 0 0 ; zeros(1,13)].*meV; 

% Find a minimum

SysOut = INS_IEP(Sys0,Sys0,Exp);

%Note that this calculated system is not the same as the one used to
%simulate teh eigenvalues
SysOut.Output.ErrorAtDeflatedPoint
%%
clear Opt
Opt.NDeflations = 2;
Opt.Verbose = true;
% Opt.Scaled = true;
Opt.Sigma = 1e-10;
% Opt.Regularisation = 1e-8;
Opt.StepTolerance = 1e-4;
% Opt.Theta = 1;
% Opt.Alpha0 = 1e-1;
Opt.C1= 0;
% Opt.Epsilon = 0.0001;
Opt.IEPType = "Classic";

SysOut = INS_IEP(Sys0,Sys0,Exp,Opt);



%%
%One option to try if the method doesn't converge is to scale the problem
%by the initial guess given:
Opt.Scaled = true;
SysOut = INS_IEP(Sys0,Sys0,Exp,Opt);
% Note that this now causes a rank deficient jacobian. We have already
% discussed the method of regularisation to deal with this, an alternative
% approach is to change the linear solver that is used:
Opt.LinearSolver = "lsqminnorm";
SysOut = INS_IEP(Sys0,Sys0,Exp,Opt);

% There are two ways to combat this problem, use the Riemannian Gradient
% descent Lift and projection method, with an increased number of iterations:
%%
clear Opt
Opt.Method = "RGD_LP";
Opt.Verbose = true;
Opt.Scaled = true;
Opt.MaxIter = 1e5;
Opt.C1=0;
Opt.Regularisation = 1e-2;
Opt.LinearSolver = "lsqminnorm";

SysOut = INS_IEP(Sys0,Sys0,Exp,Opt)
SysOut.Output

%%
% We have already discussed using regularisation to deal with a rank 
% deficient jacobian, an alternative approach is to change the linear
% solver that is used to calculate the next step:
clear Opt
% Opt.Verbose = true;
Opt.LinearSolver = "lsqminnorm";
% Opt.C1 = 1e-18;
Opt.Scaled = true;
Opt.NDeflations = 100;
% Opt.IEPType = "Classic";
SysOut = INS_IEP(Sys0,Sys0,Exp,Opt);
% SysOut.Output
% Opt.StepTolerance = 1e-6;



