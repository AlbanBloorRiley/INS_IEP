%%CrN Example
clear SysProblem SysModel Exp
%Simulate experimental data
rcm = 29979.2458;   
meV = rcm*8.065;  
N_electrons = 10;    %Number of Cr atoms in chain
N_eigenvalues = 21; 
% rng(2)
%% Set up Example Problem 

S=1/2;  SysProblem.S = [S];    %set up Spins
Jval =1.46; SysProblem.J = []; 
magnitude = 1;%set up exchange terms
for i = 2:N_electrons   %Loop over all electrons
    SysProblem.S = [SysProblem.S,S];
    SysProblem.J = [(Jval*10^((rand-0.5)*2*magnitude))*meV,zeros(1,i-2),SysProblem.J];
end

%Calculate simulated eigenvalues
[~,EE] = eigs(ham(SysProblem,[0,0,0],'sparse'),N_eigenvalues+10,'smallestreal');
EE=diag(EE);
Exp.ev=EE(1:N_eigenvalues) - EE(1);
% SysProblem.J


%% Set up model to be optimised
magnitude = 0.5;
SysModel.S = SysProblem.S;    %set up Spins
SysModel.J = SysProblem.J.*(10.^((rand(1,length(SysProblem.J))-0.5)*2*magnitude));    %set up exchange terms

Vary = SysModel; %This will vary all non-zero parameters
%% Solve problem using IEP method

%Optimse over parameters given
Opt = struct('NMinima',2,'Method','Gauss-NewtonT1','Linesearch','Basic',...
    'MaxIter',600,'theta',2,'StepTolerance',1e-6,'GradientTolerance',1e-1,...
    'ObjectiveTolerance',1e-1,'Minalpha',1e-10,'Scaled',1,...
    'deflatelinesearch',1,'IEPType','Difference','Verbose',1,'tau',0.5);

[SysOut, NIter, Flags, Iters, FinalError]= INS_IEP(SysModel,Vary,Exp,Opt);
SysModel.J

FinalError

% SysOut.J
% SysProblem.J
%% Solve problem using generic NLLS


% (startparameters,LB,UB,step)
[constants.A,constants.A0,x0,Ops,SysFixed] = Sys_Input(SysModel,Vary);
constants.ev = Exp.ev; constants.ED = 'eigs'; IEPType = "Difference";
NLLSR = @(x)NLLSresidual(x,constants,IEPType);
% LB = 0.01*x0;
% UB = 100*x0;
LB = [];UB = [];
smalleststep = 0.01;
fittype = 'single';
fittype = 'multi';
numstarts = 3;
problem = createOptimProblem('lsqnonlin','x0',x0,'objective',NLLSR,'lb',LB,'ub',UB);
problem.options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','DiffMinChange',smalleststep);
clear allmins
tic
switch fittype
    case 'single'
        [xmulti,f,output,flag] = lsqnonlin(problem);
        allmins.X = xmulti;
    case 'multi'
        ms = MultiStart('PlotFcns',@gsplotbestf,'StartPointsToRun','bounds','Display','iter');
        [xmulti,f,flag,~,allmins] = run(ms,problem,numstarts);
end
toc
f


% %%
% %Set up Mint variables - only works for Nelectrons = 6.
% % Mint still runs with less than 6 electrons? Just ignores unncessary
% % coordinates?
% a = [0;-3.3807;-5.7290;-5.5960;-3.1148;0.2482];
% b = [0;0.0840;-2.3170;-5.6470;-7.8270;-7.7410];
% c = [0; -0.2765 ; 0.1103 ; -0.2956 ; 0.2920 ; -0.1803];
% SysMint.Coords = [a b c];
% SysMint.FormFactor = 'Cr3';
% MintOpt.NumEigs = 24;
% ExpMint.SpectrumType = 'SE';
% ExpMint.Temperature = [1.5, 6, 15];
% ExpMint.Energy = -1:0.001:4;
% ExpMint.lwfwhm = 0.2;
% ExpMint.Q = [0.1,0.5,1,1.5,2]; 
% %have to 'press any key' in the command window to move between minima
% for i = 1:length(SysOut)
%     SysMint = SysOut(i);
%     SysMint.Coords = [a b c]; SysMint.FormFactor = 'Cr3';
%     [cross,Eigs] = mint(SysMint,ExpMint,MintOpt);
%     plot(ExpMint.Energy,cross)
%     xlim([0 3.25])
%     xlabel('Energy (meV)')
%     ylabel('Intensity (arb. units)')
%     legend('1.5 K','6 K','15 K')
%     if i <length(SysOut)
%         pause
%     end
% end







%% 
function r = NLLSresidual(x,constants,IEPType)
if IEPType == "Classic"
    [~,r] = IEP_Evaluate_full(x,constants);
elseif IEPType == "Difference"
    [~,r] = IEP_Evaluate_diff(x,constants);
end
end


% function to simulate the experimental eigenvalues
function [Exp,Sys] = Cr_solution(N_electrons,Neigenvalues)
rcm = 29979.2458;   meV = rcm*8.065;    %conversions
% B20 = (-0.041/3).*meV; B22 = 0.007.*meV; % Known stevens parameters
S=3/2;  Sys.S = [S];    %set up Spins
% B2 = [B22 0 B20 0 0];   Sys.B2 = [B2];  %set up Stevens parameters
Jval =1.46*meV; Sys.J = [];    %set up exchange terms
for i = 2:N_electrons   %Loop over all electrons
    Sys.S = [Sys.S,S];
%     Sys.B2 = [Sys.B2;B2];
    Sys.J = [Jval,zeros(1,i-2),Sys.J];
end
Sys.J = [ Jval/2 0 0 0 0 Jval 0 0 0 Jval 0 0 Jval 0 Jval/2];


%Calculate simulated eigenvalues
[~,EE] = eigs(ham(Sys,[0,0,0],'sparse'),Neigenvalues+10,'smallestreal');
EE=diag(EE);
Exp.ev=EE(1:Neigenvalues) - EE(1);
end




%%%


%uniform random number within order of magnitude
%take out B22 B20, just have J values, Vary spin centres. 
% keep jvalues within 1 order of magnitude in meV.
%isotropic then anisotropic
%could do 12 1/2 spins. 

