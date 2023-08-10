%Define unit conversions
rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
Tesla = meV*0.116;   % Tesla         to MHz

%% Mn12
clear Sys
Sys1 = Mn12_Spin_Sys_3(1,1);
Sys.S=Sys1.S;
Sys.B2 = [100,0,-1000,0,0];
Sys.B4 = [-1,0,0,0,-1,0,0,0,0];
Exp.ev=Sys1.ev;
Vary=Sys;
%% Cr_n
clear Sys
Sys1=Cr_Spin_Sys_3(6);
r=1;
% Sys.B2 = [1;1;1;1]*[1000 0 -1000 0 0];
Sys.S=Sys1.S;
Sys.B2=round(Sys1.B2,r,'significant');
Sys.ee=round(Sys1.ee,r,'significant');
Exp.ev=Sys1.ev;
Vary = Sys;
%% Test 1
clear Sys Sys1
rng(0)
NumEigs = 10;
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
rng(0)
[Sys,Vary,Exp]=Dy_Sys;
%%
clear Sys Vary
[Sys,Vary,Exp]=Mn6_Sys;

Opt = struct('NMinima',1,'Method','Gauss-NewtonT1','Linesearch','Armijo',...
    'MaxIter',100,'p',2,'StepTolerance',1e-4,'GradientTolerance',1e-1,...
    'ObjectiveTolerance',1e-1,'Minalpha',1e-18,'Scaled',1,...
    'deflatelinesearch',1,'IEPType','Classic','Verbose',1,'tau',0.5);
%Now find the minima
[SysOut, NIter, Flags, Iters, FinalError]= INS_IEP(Sys,Vary,Exp,Opt)
%%






Opt = struct('NMinima',1,'Method','Gauss-NewtonT1','Linesearch','Quadratic',...
    'MaxIter',1000,'p','exp','StepTolerance',1e-6,'GradientTolerance',1e-1,...
    'ObjectiveTolerance',1e-1,'Minalpha',1e-18,'Scaled',1,...
    'deflatelinesearch',1,'IEPType','Classic','Verbose',0,'tau',0.5);
[SysOut1, NIter1, Flags, Iters, FinalError]= INS_IEP(Sys,Vary,Exp,Opt)




%% 

clear SysOut NIter  
j=1;    idx=[];
for i = 1:length(NIter1)
    if Flags{i} == "Objective less than tolerance"||Flags{i} == "Gradient less than tolerance"
            idx=[idx,i];
%         SysOut(j) =SysOut1(i);
%         NIter(j)=NIter1(i);
%         j=j+1;
    end
end
SysOut=SysOut1(idx)
NIter=NIter1(idx);
FinalError(idx)




