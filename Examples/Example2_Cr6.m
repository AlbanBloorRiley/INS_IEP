%% Cr_6 Chromium-6 Chains [1]
clear all
rcm = 29979.2458;   meV = rcm*8.065;  %some conversions
%% Set up Hamiltonian model
N_electrons = 6; %Can vary the length of the chain
%Spins of 3/2
S=3/2;   Sys.S = ones(1,N_electrons)*S; 
% Use two stevens operators (D and E)
Sys.B2 = [ones(N_electrons,1) zeros(N_electrons,1) ones(N_electrons,1)...
    zeros(N_electrons,1) zeros(N_electrons,1) ];
%Because there is more than one spin centre we include the 
% isotropic electron-electron exchange coupling terms for nearest 
% neighbours using by using Sys.J (can also use Sys.ee).
J = []; %ee= [];
for i = 2:N_electrons
    J = [1,zeros(1,i-2),J];
    % ee = [eye(3);zeros(3*(i-2),3);ee];
end
Sys.J = J;
% Sys.ee = ee;

%% Simulate the eigenvalues:
ExpSys.S = Sys.S;
ExpSys.J = Sys.J*1.46*meV;
ExpSys.B2 = Sys.B2.*[0.007.*meV 0 (-0.041/3).*meV 0 0];
H = ham(ExpSys,[0,0,0]);
[~,E]=eig(H); EE = diag(E);  Exp.ev=EE(1:24)-EE(1);

Exp.ev = [0,0.355,0.457,0.497,1.576,1.577,1.592,1.629,1.632,2.97,2.98,3.002,3.004,3.01,3.038,3.821,3.824,3.827,3.837,3.856,3.879,3.888,3.895,3.903];
%% Set up initial guess 
% Since it is known a priori that each spin centre will have the same value parameters we will pin the parameters here, by setting the initial gues as the same value
Sys1.S = Sys.S;
Sys1.B2 = Sys.B2.*[1 0 -1 0 0];
% Sys1.B2 =  [1     0    -1     0     0;
%            1     0    -1     0     0;
%            1     0    -1     0     0;
%            1     0    -1     0     0;
%            1     0    -1     0     0;
%            1     0    -1     0     0];
Sys1.J = Sys.J*1e2;
% Sys1.J = [100,0,0,0,0,100,0,0,0,100,0,0,100,0,100];


Vary1 = Sys;

% Calculate the minimising system: 
SysOut1= INS_IEP(Sys1,Vary1,Exp);


%% Using Sys.ee
% Can alternatively Sys.ee to describe the total coupling matrix, this can
% be Nx3 or 3Nx3 if the antisymmetric terms are also required. 
Sys2 = Sys1;
ee = [];
for i = 2:N_electrons
    ee = [eye(3);zeros(3*(i-2),3);ee];
end
Sys2.ee = ee.*1e2;

%Note that it is not possible to use both  Sys.J and Sys.ee at the same time
%When Sys.ee is Nx3 parameters can be pinned as before by inputting the
%same initial guess values, when Sys.ee is 3Nx3 this can be done for the
%diagonal entries and off diagonal entries if the input guess is skew
%symmetric. 

Sys2 = rmfield(Sys2,'J');
%and Vary has to be updated: 
Vary2 =Sys2;

SysOut2= INS_IEP(Sys2,Vary2,Exp);

%Check that the output structures are equivalent:
all(all(ham(SysOut1,[0,0,0]) == ham(SysOut2,[0,0,0])))
%%
%It is also possible to find multiple numerical solutions to this problem,
%even if they make no physical sense. In this case it is necessary to
%change one of the deflation parmeters to make the deflation a little more 
% sensitive. 

Opt = struct('NDeflations',4,'Sigma',1e-7);
SysOut= INS_IEP(Sys1,Vary1,Exp,Opt);



% [1] - Michael L. Baker et al. “Varying spin state composition by the choice of capping
% ligand in a family of molecular chains: detailed analysis of magnetic properties
% of chromium(iii) horseshoes”. en. In: Dalton Transactions 40.12 (2011), p. 2725.
% issn: 1477-9226, 1477-9234. doi: 10.1039/c0dt01243b. url: http://xlink.
% rsc.org/?DOI=c0dt01243b.