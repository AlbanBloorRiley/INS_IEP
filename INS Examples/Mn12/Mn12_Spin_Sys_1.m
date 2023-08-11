function Mn12_Spin_Sys
%conversions
rcm = 29979.2458;    % reciprocal cm to MHz
meV = rcm*8.065;     % meV           to MHz

B20 = -0.0570*meV/3; %(D = 3*B02)
B40 = (-2.78*10^-6)*meV;
B44 = (-3.2*10^-6)*meV; 
B22 = (6.8*10^-4)*meV; %(=E)

S = 10; 
Sys.S = S; 

n = (2*S+1)^1;
Sys.B2 = [B22 0 B20 0 0];        % B(k=2,q) with q = +2,+1,0,-1,-2
Sys.B4 = [B44 0 0 0 B40 0 0 0 0];  % B(k=4,q) with q = +4,+3,+2,+1,0,-1,-2,-3,-4
H = sham(Sys,[0,0,0]); 
[Vecs,E]=eig(H);
EE = diag(E);




    A0 = sparse(n,n);
    %The eigenvalue Shift
    A{1} = speye(n);    
    %The Stevens Operators 
    A{2} = stev(S,2,0);
    A{3} = stev(S,4,0);
    A{4} = stev(S,4,4);
    A{5} = stev(S,2,2);
assignin('base','EE',EE)
assignin('base','A',A)
assignin('base','A0',A0)

%%
Exp.SpectrumType = 'SE';
Exp.Temperature = 24;
Exp.Energy = -1.5:0.005:1.5; Exp.lwfwhm = 0.1; Exp.Q = 0.2:0.05:1;
Sys.FormFactor = 'Mn3';
Sys.Coords = [0 0 0];
[cross,Eigs] = mint(Sys,Exp);

figure
plot(Exp.Energy,cross)