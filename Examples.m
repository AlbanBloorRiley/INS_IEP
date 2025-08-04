%% Mn12 Manganese-12-acetate [1,2,3]
clear all
rcm = 29979.2458;   meV = rcm*8.065;  %some conversions

%% Simulate the eigenvalues based on results from [1]
B20 = -0.0570*meV/3; B40 = (-2.78*10^-6)*meV;
B44 = (-3.2*10^-6)*meV;  B22 = (6.8*10^-4)*meV; %(=E)
% B44 = (5.1*10^-6)*meV;   B22 = (5.3*10^-4)*meV;  %the other solution
ExpSys.S = 10;   ExpSys.B2 = [B22 0 B20 0 0];        % B(k=2,q) with q = +2,+1,0,-1,-2
ExpSys.B4 = [B44 0 0 0 B40 0 0 0 0];  % B(k=4,q) with q = +4,+3,+2,+1,0,-1,-2,-3,-4
H = ham(ExpSys,[0,0,0]);
% Exp.ev = [0,1.46e-13,1.24,1.24,2.3,2.3,3.18,3.18,3.91,3.91,4.5,4.5,4.97,4.97,5.32,5.32,5.54,5.59,5.69,5.75,5.78].*meV;
[~,E]=eig(H); EE = diag(E);  Exp.ev=EE-EE(1);

%% Set up the model spin system and initial values

Sys0.S=10;  %Only one spin centre is used(Giant Spin approximation)

%Use 4 Stevens operators:
Sys0.B2 = [-100,0,-1000,0,0];   
Sys0.B4 = [-1,0,0,0,-1,0,0,0,0];

%Vary all non zero parameters (no Fixed parameters):
Vary=Sys0; 

%Calculate a minimising system using INS_IEP:
SysOut1= INS_IEP(Sys0,Vary,Exp);

%The value of the objective function can be recoved
SysOut1.Output.ErrorAtDeflatedPoint
%Clearly a solution has been found. 

%% Alternate IEP formulation
% Note that there are two formulations of the Inverse eigenvalue problem,
% this is due to the fact that INS experiments caluculate the difference,
% between eigenvalues not their explicit values. Thus the "Difference" type
% matches to this difference directly (reducing the number of equations
% used to match the data by one), whearas the "Classic" type adds an
% additional parameter to fit (the value of the groundstate/smallest
% eigenvalue. In most cases both formulations will find a minimum:

%this alternate formulation is selected by using the optional input
%structure Opt:

Opt.IEPType = "Classic";
SysOut2= INS_IEP(Sys0,Vary,Exp,Opt);

%Note that his method actually finds a different solution than the one 
% above it is reflected in the plane B22 = 0;

%% Deflation and multiple solutions
%As can be seen in figure 1 this IEP has 4 distinct solutions. Deflation, 
% [4] is a technique to systematically find multiple local minima of a system. 

%The syntax is simple:
Opt.NDeflations = 3;
SysOut3= INS_IEP(Sys0,Vary,Exp,Opt);
%This finds three different solutions to the system

%% Known solutions
%If there are other minima to the system that are already know or have been
%previously calculated they can also be input:
  
Opt = struct('NDeflations',5,'SysFound',SysOut3);
SysOut4= INS_IEP(Sys0,Vary,Exp,Opt);
%This finds all 4 solutions to the system, using the three previously found
% inima. 



%% Cr_6 Chromium-6 Chains [5]
clear all
rcm = 29979.2458;   meV = rcm*8.065; 
%set up model
N_electrons = 6; %Original paper has 6 chains, but can work for other lengths
%Spins of 3/2
S=3/2;   Sys.S = ones(1,N_electrons)*S; 
% Use two stevens operators (D and E)
Sys.B2 = [ones(N_electrons,1) zeros(N_electrons,1) ones(N_electrons,1)...
    zeros(N_electrons,1) zeros(N_electrons,1) ];
%Include isotropic electron-electron exchange coupling terms for nearest 
% neighbours using Sys.J
J = []; %ee= [];
for i = 2:N_electrons
    J = [1,zeros(1,i-2),J];
    % ee = [eye(3);zeros(3*(i-2),3);ee];
end
Sys.J = J;
% Sys.ee = ee;
%Simulate the eigenvalues:
ExpSys.S = Sys.S;
ExpSys.J = Sys.J*1.46*meV;
ExpSys.B2 = Sys.B2.*[0.007.*meV 0 (-0.041/3).*meV 0 0];
H = ham(ExpSys,[0,0,0]);
[~,E]=eig(H); EE = diag(E);  Exp.ev=EE(1:24)-EE(1);


%%
%Set up initial guess, using Sys.J
Sys1.S = Sys.S;
Sys1.B2 = Sys.B2.*[1 0 -1 0 0];
Sys1.J = Sys.J*1e2;
Vary1 = Sys;

SysOut1= INS_IEP(Sys1,Vary1,Exp);

%Can alternatively Sys.ee to describe the total coupling matrix, if we only
%include the isotropic terms then the end result will be the same:
Sys2 = Sys1;
ee = [];
for i = 2:N_electrons
    ee = [eye(3);zeros(3*(i-2),3);ee];
end
Sys2.ee = ee.*1e2;

%Note that it is not possible to use both  Sys.J and Sys.ee at the same time
Sys2 = rmfield(Sys2,'J');
%and Vary has to be updated: 
Vary2 =Sys2;
SysOut2= INS_IEP(Sys2,Vary2,Exp);

%Check that the output structures are equivalent:
all(all(ham(SysOut1,[0,0,0]) == ham(SysOut2,[0,0,0])))

%It is also possible to find multiple numerical solutions to this problem,
%even if they make no physical sense. In this case it is necessary to
%change one of the deflation parmeters. 
Opt = struct('NDeflations',4,'Sigma',1e-7);
SysOut= INS_IEP(Sys1,Vary1,Exp,Opt);





%% Dy_6 Dysprosium
clear all
rcm = 29979.2458;   meV = rcm*8.065; 
%Experimental eigenvalues
Exp.ev = [0;0;0.72;0.72;1.12;1.12;1.28;1.28].*meV;

Sys.S = [1/2 15/2 1/2];
% B20 = -0.173585036568371*meV;B22 = 0.443886669455654*meV;
% B40 = -0.000767621238475*meV;B60 = -0.000125378265669*meV;
% Jex1 = 0.076588615657158*meV;Jex2 = 0.048947841762338*meV;
% 
% B20 = -0.173585*meV;B22 =  0.443886*meV;
% B40 = -0.000767*meV;B60 = -0.000125*meV;
% Jex1 = 0.076588*meV;Jex2 = 0.048947*meV;

B20 = -0.1*meV;B22 =  0.1*meV;
B40 = -0.001*meV;B60 = -0.0001*meV;
Jex1 = 0.1*meV;Jex2 = 0.01*meV;
Sys.B2 = [0 0 0 0 0 ; B22 0 B20 0 0 ; 0 0 0 0 0 ];
% Sys.B2 = [0 0 0 0 0 ; 0 0 B20 0 0 ; 0 0 0 0 0 ];
Sys.B4 = [0 0 0 0 0 0 0 0 0; 0 0 0 0 B40 0 0 0 0; 0 0 0 0 0 0 0 0 0];
Sys.B6 = [0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 B60 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0];
% Sys.J = [Jex1 0 Jex2];

% Sys.ee = [diag([Jex1,Jex1,Jex1]);zeros(3);diag([Jex2,Jex2,Jex2])];
% 
Sys.ee = [diag([Jex1,Jex1,Jex1+1]);zeros(3);diag([Jex2,Jex2,Jex2+1])];
% Vary.ee = [[1,0,1;0,1,0;1,0,1];zeros(3);[1,0,1;0,1,0;1,0,1]];
Vary = Sys;
%
% Vary.B2([2,8]) = 0;
% Vary.B4(14) = 0;
% Vary.B6(20) = 0;

% This model of the system is rank deficien and does not converge within 
% 1000 iterations with default settings:
SysOut1 = INS_IEP(Sys,Vary,Exp);
%Note the warning about rank deficiency, now if we include a regularisation
%term the method converges very quickly:
%
% Opt.Verbose = true;
Opt.Regularisation = 1e-4;
SysOut2 = INS_IEP(Sys,Vary,Exp,Opt);

% Antiisotropic? off diagonals

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



%% Example 2
% PHYSICAL REVIEW B 88, 064405 (2013)
clear all
% some constants
rcm = 29979.2458; % reciprocal cm to MHz
meV = rcm*8.065; % meV to MHz
% All parameters come from PHYSICAL REVIEW B 88, 064405 (2013)
%Setup of J=6 for Tb3+ and S=1/2 for Cu
ExpSys.S = [6 1/2];

% Exchange coupling
Jxx = -0.54;
Jyy = -0.54;
Jzz = -0.298;
ExpSys.ee = [Jxx Jyy Jzz].*meV;
%Stevens parameters
B20 = -0.197;
B22 = 0.066;
B40 = 1.5e-4;
B44 = 5e-4;
B60 = -1.25e-5;
B64 = 1.17e-4;
ExpSys.B2 = [B22 0 B20 0 0 ; zeros(1,5)].*meV;
ExpSys.B4 = [B44 0 0 0 B40 0 0 0 0 ; zeros(1,9)].*meV;
ExpSys.B6 = [0 0 B64 0 0 0 B60 0 0 0 0 0 0 ; zeros(1,13)].*meV; 
H = ham(ExpSys,[0,0,0]);
EE = eig(H);
Exp.ev = EE(1:10)-EE(1);

Sys0.S = [6 1/2];
Sys0.ee = [-0.5 -0.5 -0.25].*meV;
Sys0.B2 = [0.1 0 -0.1 0 0 ; zeros(1,5)].*meV;
Sys0.B4 = [1e-4 0 0 0 1e-4 0 0 0 0 ; zeros(1,9)].*meV;
Sys0.B6 = [0 0 1e-4 0 0 0 -1e-4 0 0 0 0 0 0 ; zeros(1,13)].*meV; 

%This is another highly nonlinear problem, with an often poorly conditioned
% Jacobian matrix, meaning the method doesnt even take a single iteration 
% as the line search terminates:
SysOut = INS_IEP(Sys0,Sys0,Exp);
SysOut.Output.ErrorAtDeflatedPoint
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
% Opt.Verbose = true;
Opt.MaxIter = 1e5;
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





%%
clear all
% some constants
rcm = 29979.2458; % reciprocal cm to MHz
meV = rcm*8.065; % meV to MHz
% All parameters come from PHYSICAL REVIEW B 88, 064405 (2013)
%Setup of J=6 for Tb3+ and S=1/2 for Cu
ExpSys.S = [6 1/2];

% Exchange coupling
Jxx = -0.54;
Jyy = -0.54;
Jzz = -0.298;
ExpSys.ee = [Jxx Jyy Jzz].*meV;
%Stevens parameters
B20 = -0.197;
B22 = 0.066;
B40 = 1.5e-4;
B44 = 5e-4;
B60 = -1.25e-5;
B64 = 1.17e-4;
ExpSys.B2 = [B22 0 B20 0 0 ; zeros(1,5)].*meV;
ExpSys.B4 = [B44 0 0 0 B40 0 0 0 0 ; zeros(1,9)].*meV;
ExpSys.B6 = [0 0 B64 0 0 0 B60 0 0 0 0 0 0 ; zeros(1,13)].*meV; 
H = ham(ExpSys,[0,0,0]);
EE = eig(H);
Exp.ev = EE(1:end)-EE(1);

%Not every parameter needs to be fitted. Only the non zero values in the
%structure 'Vary' will be varied. All the other parameters will stay fixed.
%In this example we fix the zero field splitting parameters and 
Sys0 = ExpSys;
Sys0.ee = [-1 -1 -2].*meV;
Vary.ee = Sys0.ee;
Vary.B2 = zeros(size(ExpSys.B2));
Vary.B4 = zeros(size(ExpSys.B4));
Vary.B6 = zeros(size(ExpSys.B6));
% Sys0.S = [6 1/2];
% Sys0.ee = [-0.5 -0.5 -0.25].*meV;
% Sys0.B2 = [0.1 0 -0.1 0 0 ; zeros(1,5)].*meV;
% Sys0.B4 = [1e-4 0 0 0 1e-4 0 0 0 0 ; zeros(1,9)].*meV;
% Sys0.B6 = [0 0 1e-4 0 0 0 -1e-4 0 0 0 0 0 0 ; zeros(1,13)].*meV; 

%This is another highly nonlinear problem which fails to converge with
% standard options :
SysOut = INS_IEP(Sys0,Vary,Exp);
SysOut.Output.ErrorAtDeflatedPoint
%%
%One option is to use a different numerical method:
Opt.Method = "RGD_LP";
SysOut = INS_IEP(Sys0,Vary,Exp,Opt)



% Another option is to try and find the minimum after multiple deflations.
% When using deflation is it often adventagous to scale the variables:
clear Opt
% Opt.IEPType = "Difference";
Opt.Scaled = true;
Opt.NDeflations = 5;
Opt.C1 = 1e-10;
SysOut = INS_IEP(Sys0,Sys0,Exp,Opt);
% In this case it is the 5th found system that we want.
SysOut(5)

%%

clear Opt
Opt.Method = "Newton";
Opt.Scaled = true;
Opt.NDeflations = 300;
Opt.Verbose = true;
% Opt.LinearSolver = "lsqminnorm";
% Opt.Regularisation = 1e-18;
Opt.FunctionTolerance = 1e-3;
Opt.ConvergenceFlags = "Gradient less than tolerance";
SysOut = INS_IEP(Sys0,Sys0,Exp,Opt);
% In this case it is the 5th found system that we want.
% SysOut.Output

%%
%% Longer Chromium Chains [5]
clear all
rcm = 29979.2458;   meV = rcm*8.065; 
%set up model
N_electrons = 8; %Original paper has 6 chains, but can work for other lengths
%Spins of 3/2
S=3/2;   Sys.S = ones(1,N_electrons)*S; 
% Use two stevens operators (D and E)
Sys.B2 = [ones(N_electrons,1) zeros(N_electrons,1) ones(N_electrons,1)...
    zeros(N_electrons,1) zeros(N_electrons,1) ];
%Include isotropic electron-electron exchange coupling terms for nearest 
% neighbours using Sys.J
J = []; %ee= [];
for i = 2:N_electrons
    J = [1,zeros(1,i-2),J];
    % ee = [eye(3);zeros(3*(i-2),3);ee];
end
Sys.J = J;
% Sys.ee = ee;
%Simulate the eigenvalues:
ExpSys.S = Sys.S;
ExpSys.J = Sys.J*1.46*meV;
ExpSys.B2 = Sys.B2.*[0.007.*meV 0 (-0.041/3).*meV 0 0];
tic
% H = ham(ExpSys,[0,0,0],'sparse');
[A,A0,x0,Ops,SysFixed] = Sys_Input(ExpSys,ExpSys);
H = FormA(x0,A,A0);
[~,E]=eigs(H,30,'smallestreal'); EE = diag(E);  Exp.ev=EE(1:24)-EE(1);
toc

%%
%Set up initial guess, using Sys.J
Sys0.S = Sys.S;
Sys0.B2 = Sys.B2.*[1 0 -1 0 0];
Sys0.J = Sys.J*1e2;
Vary1 = Sys;
Opt.Verbose = true;
tic
SysOut= INS_IEP(Sys0,Vary1,Exp,Opt);
toc

%%
clear Sys0
rcm = 29979.2458;   meV = rcm*8.065; 
SysSol.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2

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

SysSol.J = [J_AB J_A1 J_A2 J_A3 J_A4 J_B1 J_B2 J_B3 J_B4 J_12 J_13 J_14 J_23 J_24 J_34].*(-2); %-2J formalism

DIII = -0.029*meV; % MnIII anisotropy. Free Value 
EIII = 0*DIII;   % MnIII rhombicity. For INS, assume 0.
B20III = 3*DIII; B22III = EIII; %Converting to Stevens Operator formalism

DII = -0.000*meV; % MnII anisotropy. For INS, assume 0. 
EII = DII*0;      % MnII rhombicity. For INS, assume 0.
B20II = 3*DII; B22II = EII; %Converting to Stevens Operator formalism

SysSol.B2 = [B22III 0 B20III 0 0;
          B22III 0 B20III 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0];

%%
Sys0.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2

J_S4_S4   = -5*meV;    % MnIII - MnIII.                       
J_S4_S5_1 = 0.6*meV;  % MnIII - MnII, MnIII JT involved.     
J_S4_S5_2 = -0.2*meV; % MnIII - MnII, MnIII JT not involved. 
J_S5_S5   = -0.01*meV;  % MnII  - MnII.       
J_AB = J_S4_S4;
J_A1 = J_S4_S5_1; J_A3 = J_S4_S5_1; J_B2 = J_S4_S5_1; J_B4 = J_S4_S5_1; %For these, the MnIII JT axis is involved. Can be FM or AFM
J_A2 = J_S4_S5_2; J_A4 = J_S4_S5_2; J_B1 = J_S4_S5_2; J_B3 = J_S4_S5_2; %For these, the MnIII JT axis is NOT involved. Can only be AFM
J_12 = J_S5_S5;   J_34 = J_S5_S5; 
J_14 = 0;         J_23 = 0;% 0.01*meV; %Assumed zero. Only interacts via a pivalate bridge.
J_13 = 0;         J_24 = 0; %Assumed zero. No interaction pathway

Sys0.J = [J_AB J_A1 J_A2 J_A3 J_A4 J_B1 J_B2 J_B3 J_B4 J_12 J_13 J_14 J_23 J_24 J_34].*(-2); %-2J formalism

DIII = -0.029*meV; % MnIII anisotropy. Free Value 
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

Vary = Sys0; %Vary.J(1)=0;


Exp.ev = [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.5218, 2.5218, 2.5218, 2.5218]'.*meV;   %or
% Exp.ev = [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.316, 2.316,  2.5218, 2.5218]'.*meV;

clear Opt
% Opt.ScalingMatrix = 1e2*eye(5);
Opt.Verbose = true;
% Opt.Scaled = true;
Opt.C1 =1e-2;
Opt.Alpha0 = 1e-0;
Opt.Linesearch = "Quadratic";
% Opt.Method = "Classic";
Opt.MaxIter = 500;
% Opt.Method = "RGD_LP";
SysOut = INS_IEP(Sys0,Vary,Exp,Opt)

%%






%% General tips
%
% If failing to converge then use: scaling, regularisation, try RGD_LP or
% another method
%
% If converging at 0/1 iterations, recuce C1
%
% If rank deficiency try: Regularisation or change the LinearSolver to
% "lsqminnorm" or custom Solver. 
%
% If deflation is not effective use: Scaling, increase Theta, decrease
% Sigma
%
% 
%
%



%% Test 1
clear Sys Sys1 Exp
rng(1)
NumEigs = 100;
Sys1.S = [2.5 2 2];
B20I = 6*rcm;   B22I = 0.1*rcm;
B20II = 15*rcm; B22II = 0.5*rcm;
J1=15*rcm; J2=5*rcm;
Sys1.B2 =[B22I 0 B20I 0 0; B22II 0 B20II 0 0; B22II 0 B20II 0 0];
Sys1.J = [J1 J1 J2];
EE=eig(ham(Sys1,[0 0 0]));
ev=EE(1:NumEigs)-EE(1);
Sys.S=Sys1.S;
B22I = (0.1+0.05*randn)*rcm;   B20I = (6+1*randn)*rcm; 
B22II = (0.5+0.1*randn)*rcm;  B20II = (15+5*randn)*rcm; 
J1=(15+1*randn)*rcm; J2=(5+1*randn)*rcm;
Sys.B2 =[B22I 0 B20I 0 0; B22II 0 B20II 0 0; B22II 0 B20II 0 0];
Sys.J = [J1 J1 J2];
Vary = Sys;
Exp.ev=ev;
Opt = struct('NDeflations',2,'Method','Good_GN','verbose',true,'Linesearch','Quadratic','Scaled',true,'c1',1e-10);
SysOut= INS_IEP(Sys,Vary,Exp,Opt)





% [1] - Roland Bircher et al. “Transverse magnetic anisotropy in Mn 12 
% acetate: Direct determination by inelastic neutron scattering”. en. In: 
% Physical Review B 70.21 (Dec. 2004), p. 212413. issn: 1098-0121, 1550-235X

% [2] - J. R. Friedman, M. P. Sarachik, J. Tejada, and R. Ziolo, 
% Macroscopic Measurement of Resonant Magnetization Tunneling in High-Spin 
% Molecules, Phys. Rev. Lett., 76 (1996), pp. 3830–3833, https://doi.org/10.
% 1103/physrevlett.76.3830, https://link.aps.org/doi/10.1103/PhysRevLett.76.3830. 
% Publisher: American Physical Society (APS).

% [3] - R. Sessoli, D. Gatteschi, A. Caneschi, and M. A. Novak, Magnetic 
% bistability in a metal-ion cluster, Nature, 365 (1993), pp. 141–143, 
% https://doi.org/10.1038/365141a0, https://www.nature.com/articles/365141a0.
% Publisher: Springer Science and Business Media LLC.




% [6] - A Bloor Riley RGDLP

% [4] - A Bloor Riley Deflation

% [5] - Michael L. Baker et al. “Varying spin state composition by the choice of capping
% ligand in a family of molecular chains: detailed analysis of magnetic properties
% of chromium(iii) horseshoes”. en. In: Dalton Transactions 40.12 (2011), p. 2725.
% issn: 1477-9226, 1477-9234. doi: 10.1039/c0dt01243b. url: http://xlink.
% rsc.org/?DOI=c0dt01243b.