rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
Tesla = meV*0.116;   % Tesla         to MHz


[Sys,Exp] = Mn12_Spin_Sys_3(1,6);



constants.ev = Exp.ev;
constants.ED = "eig";

% constants.ev = Exp.ev(1:5);

Vary.B2 = [1 0 0 0 0];
Vary.B4 = [1 0 0 0 0 0 0 0 0];


[AA,A0,xint,Ops,SysFixed] =Sys_Input(Sys,Vary);
A{1} = AA{1};
A{2} =AA{2};
constants.A = A;
constants.A = AA;
constants.A0 = A0;
constants.eigenvalueSD = ones(length(Exp.ev),1);
constants.eigenvalueDifferenceSD = ones(length(Exp.ev)-1,1);

x = linspace(0,2.5,500);
y = linspace(-1.5,1.5,500);
clear X Y Z
[X,Y] = meshgrid(x,y);
% X = X;
% Y = Y;

for i = 1:length(X)
    for j = 1:length(Y)
            Z(i,j) = IEP_Evaluate_full([X(i,j)*1e-3.*meV;Y(i,j)*1e-5.*meV],constants);
    end
end
%%
[C,P] = contourf(X,Y,log(Z),100);
P.LineColor = 'none';
% colormap("gray")



