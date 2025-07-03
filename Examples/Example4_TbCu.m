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
ExpSys.ee = [-0.54 -0.54 -0.298].*meV;
%Stevens parameters
ExpSys.B2 = [0.066 0 -0.197 0 0 ; zeros(1,5)].*meV;
ExpSys.B4 = [5e-4 0 0 0 1.5e-4 0 0 0 0 ; zeros(1,9)].*meV;
ExpSys.B6 = [0 0 B64 0 0 0 -1.25e-5 0 0 0 0 0 0 ; zeros(1,13)].*meV; 
H = ham(ExpSys,[0,0,0]);
EE = eig(H);
%only save a subset of the eigenvalues
% Exp.ev = EE(1:10)-EE(1);
Exp.ev = EE(1:end)-EE(1);

%% Set up the model spin system and initial values

Sys0.S = [6 1/2];
Sys0.ee = [-0.5 -0.5 -0.05].*meV;
Sys0.B2 = [0.1 0 -0.1 0 0 ; zeros(1,5)].*meV;
Sys0.B4 = [1e-10 0 0 0 1e-10 0 0 0 0 ; zeros(1,9)].*meV;
Sys0.B6 = [0 0 1e-10 0 0 0 -1e-10 0 0 0 0 0 0 ; zeros(1,13)].*meV; 

%% Find a minimum

SysOut = INS_IEP(Sys0,Sys0,Exp);

SysOut.Output.ErrorAtDeflatedPoint
