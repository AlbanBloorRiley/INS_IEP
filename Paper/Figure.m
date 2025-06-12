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

Vary = ExpSys; Vary.B2(3)=0; Vary.B4(5)=0;
[A,A0,x0,Ops,SysFixed] = Sys_Input(ExpSys,Vary);
constants.A = A;
constants.A0 = A0;
constants.ev=Exp.ev;
N=100;
x = linspace(-2.5e-3,2.5e-3,N);
y = linspace(-1.5e-5,1.5e-5,N);
[X,Y] =meshgrid(x,y);
Z=zeros(length(x),length(y));   
for i=1:length(x)
    for j=1:length(y)
        Z(i,j) = IEP_Evaluate_diff([x(j);y(i)].*meV,constants);
    end
end
%
contourf(X,Y,Z,100,'edgecolor','none')