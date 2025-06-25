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

Opt = struct('NDeflations',4,'SysFound',SysOut3);
SysOut4= INS_IEP(Sys0,Vary,Exp,Opt);
%This finds all 4 solutions to the system, using the three previously found
% inima. 


%%
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

% [4] - deflation
%%
Sys0.B2 = [-1,0,-4000,0,0];   
Sys0.B4 = [1e-2,0,0,0,-1,0,0,0,0];


Opt = struct('NDeflations',10,'Method','RGD_LP');
Opt.Scaled = true;
Opt.Sigma = 1e-8;
Opt.Theta = 5;
SysOut4= INS_IEP(Sys0,Vary,Exp,Opt);
SysOut4.B2