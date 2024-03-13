%Define unit conversions
rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
Tesla = meV*0.116;   % Tesla         to MHz

%% Mn12
clear Sys Exp
[Sys1,Exp] = Mn12_Spin_Sys_3(1,1);
Sys.S=Sys1.S;
Sys.B2 = [100,0,-1000,0,0];
Sys.B4 = [-1,0,0,0,-1,0,0,0,0];
Vary=Sys;
%% Cr_n
clear Sys Exp
[Sys1,Exp]=Cr_Spin_Sys_3(4);
H=ham(Sys1,[0,0,0],'sparse'); [Vecs,EE] = eig(full(H),'vector');
EE=EE(1:24);
Exp.ev=EE-EE(1);

r=1;
% Sys.B2 = [1;1;1;1]*[1000 0 -1000 0 0];
% Sys.S=Sys1.S;
% Sys.B2=round(Sys1.B2,r,'significant');
% Sys.ee=round(Sys1.ee,r,'significant');
% Vary = Sys;
N_electrons = 4;
B20 = (-0.041/3).*meV; % convert from meV to MHz
B22 = 0.007.*meV; % convert from meV to MHz
B20 = -1000;
B22=1000;
BValues = [-3304.4;1692.5;353000];
S=3/2;
Sys.S = [S];
for i = 2:N_electrons
    Sys.S = [Sys.S,S];
end
B2 = [B22 0 B20 0 0]; % B(k=2,q) with q = +2,+1,0,-1,-2\% Sys.B2=round(Sys.B2,r,'significant');
B2=round(B2,r,'significant');
Sys.B2 = [B2];
for i = 2:N_electrons
    Sys.B2 = [Sys.B2;B2];
end
J = 1.46*meV;   ee = [];ee1=[];JJ = [];
% J = 0.1*meV;
J = round(J,r,'significant');
for i = 2:N_electrons
%     JJ = [1+i*1e-4,zeros(1,i-2),JJ];
    JJ = [J,zeros(1,i-2),JJ];

    ee = [eye(3);zeros(3*(i-2),3);ee];
end
Sys.J = JJ;
Vary = Sys;
% Sys.ee = J.*ee;
% J=1;



Opt = struct('NMinima',3,'Method','Gauss-NewtonT1','Linesearch','Basic',...
    'MaxIter',1000,'theta',2,'StepTolerance',1e-6,'GradientTolerance',1e-1,...
    'ObjectiveTolerance',1e-1,'Minalpha',1e-10,'Scaled',1,'epsilon',0.0001,...
    'deflatelinesearch',1,'IEPType','Difference','Verbose',0,'tau',0.5);
[SysOut, NIter1, Flags, Iters, FinalError]= INS_IEP(Sys,Vary,Exp,Opt)

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
rng(0)
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
%%
clear Sys
[Sys,Vary,Exp]=Dy_Sys;
%%
clear Sys Vary
[Sys,Vary,Exp]=Mn6_Sys;

%%
Opt = struct('Eigensolver','eig','NDeflations',4,'Method','Newton','Linesearch','No',...
    'MaxIter',10,'theta',2,'StepTolerance',1e-6,'GradientTolerance',1e-1,...
    'ObjectiveTolerance',1e-1,'Minalpha',1e-10,'Scaled',false,'epsilon',0.001,...
    'deflatelinesearch',true,'IEPType','Difference','Verbose',true,'tau',0.5,"SysFound",SysFound);
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




