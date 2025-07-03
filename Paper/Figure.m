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
        [F,R,J] = IEP_Evaluate_diff([x(j);y(i)].*meV,constants);
        Z(i,j) = sqrt(F);
        % Z(i,j) = norm(J'*R);
    end
end
%%
f = figure(1);
% subplot(1,2,1)
contourf(X,Y,Z,100,'edgecolor','none');
colorbar
%
ylabel("B^4_4 meV")
xlabel("B^2_2 meV")

f.Units = 'centimeters';
f.Position = [-50 10 20 14];
%

print(f, 'Mn12_contour.png', '-dpng')

%%
Opt.NDeflations = 4;
Sys0 = ExpSys;
Sys0.B2 = [-100,0,-1000,0,0];   
Sys0.B4 = [-1,0,0,0,-1,0,0,0,0];
% Opt.Verbose = true;
SysOut = INS_IEP(Sys0, Sys0,Exp,Opt)

%%

subplot(1,2,2)
%%
f = figure(2);
hold on
cla
lgnd = string;
for i = 1:4
     x = 0:SysOut(i).Output.NIter-1;
        xx = sum((SysOut(i).Output.Iterates(:,1:end-1)-SysOut(i).Output.Iterates(:,2:end)).^2);
        semilogy(x,xx,'linewidth',1.3)
        entry = ['Deflation ', num2str(i-1)];
        lgnd = [lgnd; entry];
end
    lgnd = lgnd(2:end,:);
    if lgnd(1) == "Deflation 0"
        lgnd(1) = "Undeflated";
    end
    legend(lgnd,'Location','s','fontsize',15)
hold off

xlabel("Iterations")
ylabel("Error")
ylim([1e-10,1e5])
set(gca, 'YScale', 'log')
set(gca,'YMinorGrid','on')

f.Units = 'centimeters';
f.Position = [-50 10 20 14];
print(f, 'Mn12_convergence.png', '-dpng')

%%

% f.Units = 'centimeters';
% f.Position = [-50 10 35 14];
% %
% print(f, 'Mn12_figure.eps', '-depsc')
% print(f, 'Mn12_figure.pdf', '-dpdf')
% print(f, 'Mn12_figure.png', '-dpng')