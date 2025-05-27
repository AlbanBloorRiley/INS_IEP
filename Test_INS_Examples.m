%% Mn12 Manganese-12-acetate [1,2,3]
clear all
rcm = 29979.2458;   meV = rcm*8.065;  
%Calculate the simulated eigenvalues
B20 = -0.0570*meV/3; B40 = (-2.78*10^-6)*meV;
B44 = (-3.2*10^-6)*meV;  B22 = (6.8*10^-4)*meV; %(=E)
% B44 = (5.1*10^-6)*meV;   B22 = (5.3*10^-4)*meV;  %one of two possible solutions
ExpSys.S = 10;   ExpSys.B2 = [B22 0 B20 0 0];        % B(k=2,q) with q = +2,+1,0,-1,-2
ExpSys.B4 = [B44 0 0 0 B40 0 0 0 0];  % B(k=4,q) with q = +4,+3,+2,+1,0,-1,-2,-3,-4
H = ham(ExpSys,[0,0,0]);
[~,E]=eig(H); EE = diag(E);  Exp.ev=EE-EE(1);

%Set up the System and intial guess
%Only one spin centre (Using Giant Spin approximation):
Sys0.S=ExpSys.S;  
%Use 4 Stevens operators:
Sys0.B2 = [-100,0,-1000,0,0];   
Sys0.B4 = [-1,0,0,0,-1,0,0,0,0];
%Vary all non zero parameters (no Fixed parameters):
Vary=Sys0; 

%Calculate a minimising system:
SysOut= INS_IEP(Sys0,Vary,Exp);


%Note that there are two formulations of the Inverse eigenvalue problem,
%this is due to the fact that INS experiments caluculate the difference,
% between eigenvalues not their explicit values. Thus the "Difference" type
% matches to this difference directly (reducing the number of equations
% used to match the data by one), whearas the "Classic" type adds an
% additional parameter to fit (the value of the groundstate/smallest
% eigenvalue. In most cases both formulations will find a minimum:
Opt = struct('IEPType', 'Classic');
SysOut= INS_IEP(Sys0,Vary,Exp,Opt);

%If multiple solutions are desired can use deflation [4]:
Opt = struct('NDeflations',3);
SysOut= INS_IEP(Sys0,Vary,Exp,Opt);

%It is also possible to use previously found solutions in the deflation
%process:
Opt = struct('NDeflations',4,'SysFound',SysOut);
SysOut= INS_IEP(Sys0,Vary,Exp,Opt);

% SysOut.B2

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
B20 = -0.173585036568371*meV;B22 = 0.443886669455654*meV;
B40 = -0.000767621238475*meV;B60 = -0.000125378265669*meV;
Jex1 = 0.076588615657158*meV;Jex2 = 0.048947841762338*meV;
% 
% B20 = -0.173585*meV;B22 =  0.443886*meV;
B40 = -0.000767*meV;B60 = -0.000125*meV;
% Jex1 = 0.076588*meV;Jex2 = 0.048947*meV;

B20 = -0.1*meV;B22 =  0.4*meV;
B40 = -0.0007*meV;B60 = -0.0001*meV;
Jex1 = 0.07*meV;Jex2 = 0.04*meV;
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
SysOut = INS_IEP(Sys,Vary,Exp,Opt);
%Note the warning about rank deficiency, now if we include a regularisation
%term the method converges very quickly:
%
Opt.Regularisation = 1e-8;
SysOut = INS_IEP(Sys,Vary,Exp,Opt);

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

%In this case the standard method converges (the line search terminates), 
% but the value of the objective function is large, as well as the norm of
% the gradient so this is not a minimum, most likely due to the highly non
% linear nature of the problem:
Opt.Verbose=true
SysOutB = INS_IEP(Sys0,Sys0,ExpB,Opt)

%In this case it may be better to use an alternaitve minimisation method, 
%for example the Riemannian Gradient Descent Lift and Projection method [6]:
Opt = struct('Method', 'RGD_LP');
SysOutB = INS_IEP(Sys0,Sys0,ExpB,Opt)

Sys0.B4 = [zeros(1,8),-1e-2].*meV;
% Opt = struct('Method','RGD_LP');
SysOutA = INS_IEP(Sys0,Sys0,ExpA,Opt)




%%
ExpSys.S = Sys.S



Vary = Sys;
% Sys.ee = J.*ee;
% J=1;


% Opt = struct('NDeflations',3,'verbose',true);
% Opt = struct('NDeflations',1,'Method','LP','Linesearch','Basic',...
%     'MaxIter',1000,'theta',2,'StepTolerance',1e-6,'GradientTolerance',0,...
% 'Minalpha',1e-10,'Scaled',true,'epsilon',0.001,...
%     'IEPType','Difference','Verbose',true,'c1',1e-6);


[A] =Sys_Input(Sys,Vary);
A{end+1} = eye(size(A{1}));
    B = zeros(length(A));
    for i = 1:(length(A))
        for j = 1:(length(A))
            B(i,j) = sum(sum(A{j}'.*A{i}));
        end
    end
Opt = struct('NDeflations',1,'Method','LP','Linesearch','Armijo',...
    'MaxIter',1e5,'theta',2,'StepTolerance',1e-6,'GradientTolerance',0,...
'Minalpha',1e-10,'Scaled',false,'epsilon',0.001,...
    'IEPType','Difference','Verbose',true,'c1',1e-6,'alpha0',100);
[SysOut]= INS_IEP(Sys,Vary,Exp,Opt)

SysOut.B2



%For many minima:
% Opt = struct('NMinima',8,'Method','Newton','Linesearch','Basic',...
%     'MaxIter',1000,'theta',2,'StepTolerance',1e-6,'GradientTolerance',1e-1,...
%     'ObjectiveTolerance',1e-1,'Minalpha',1e-10,'Scaled',1,'epsilon',0.001,...
%     'deflatelinesearch',0,'IEPType','Difference','Verbose',1,'tau',0.5);
% Sys.B2=round(Sys.B2,r,'significant');
% Sys.ee=round(Sys.ee,r,'significant');
% Sys.J=round(Sys.J,r,'significant');

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











%%
clear Sys
[Sys,Vary,Exp]=Dy_Sys;
warning('Off','MATLAB:nearlySingularMatrix')
Opt = struct('Eigensolver','eig','NDeflations',10,'Method','RGD_LP','Linesearch','Armijo',...
    'MaxIter',1000,'theta',2,'StepTolerance',1e-6,'GradientTolerance',1e-10,...
    'ObjectiveTolerance',1e-6,'Minalpha',1e-18,'Scaled',true,'epsilon',0.1,...
    'IEPType','Classic','Verbose',true,'c1',1e-10);
[SysOut]= INS_IEP(Sys,Vary,Exp,Opt)
warning('On','MATLAB:nearlySingularMatrix')
%%
clear Sys Vary
[Sys,Vary,Exp]=Mn6_Sys;
%%
Opt = struct('Eigensolver','eig','NDeflations',10,'Method','Good_GN','Linesearch','Armijo',...
    'MaxIter',1000,'theta',2,'StepTolerance',1e-6,'GradientTolerance',0,...
    'ObjectiveTolerance',0,'Minalpha',1e-13,'Scaled',true,'epsilon',0.1,...
    'IEPType','Difference','Verbose',false,'tau',0.5,'c1',1e-10,...
    "MuLinesearch","No","DeflatedLinesearch","Armijo");
[SysOut]= INS_IEP(Sys,Vary,Exp,Opt)

% SysOut.B2

% SysOut.ee
%% Test scale invariance
Opt = struct('NMinima',1,'Method','Newton','Linesearch','No',...
    'MaxIter',1000,'theta',2,'StepTolerance',1e-6,'GradientTolerance',1e-1,...
    'ObjectiveTolerance',1e-1,'Minalpha',1e-10,'epsilon',0.001,...
    'deflatelinesearch',1,'IEPType','Classic','Verbose',1,'tau',0.5);
Opt.Scaled = 1;
[SysOut1, NIter1, Flags1, Iters1, FinalError1]= INS_IEP(Sys,Vary,Exp,Opt);
Opt.Scaled = 0;
[SysOut0, NIter0, Flags0, Iters0, FinalError0]= INS_IEP(Sys,Vary,Exp,Opt);
Iters1{1}.*Iters0{1}(:,1)-Iters0{1}
%% 

j=1;    idx=[];
for i = 1:length(NIter1)
    if Flags{i} == "Objective less than tolerance"||Flags{i} == "Gradient less than tolerance"
            idx=[idx,i];
%         SysOut(j) =SysOut1(i);
%         NIter(j)=NIter1(i);
%         j=j+1;
    end
end
% SysOut=SysOut1(idx)
NIter=NIter1(idx);
FinalError(idx)








%%
function ferrs = format_errs(errs)
    ferrs = string([]);
    for k = 1:length(errs)
        ferrs(k) = sprintf('%0.1e', errs(k));
    end
end





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