function Cr_Spin_Sys_1 (N_electrons)
%N_electrons = 3;clear Sys Sys1
% conversions:
rcm = 29979.2458;    % reciprocal cm to MHz
meV = rcm*8.065;     % meV           to MHz

% Setup of the Hamiltonian in the usual easyspin way
B20 = (-0.041/3).*meV; % convert from meV to MHz
B22 = 0.007.*meV; % convert from meV to MHz
BValues = [-3304.4;1692.5;353000];
% We have 6 S=3/2 spins that are coupled together
S=3/2;
Sys.S = [S];
for i = 2:N_electrons
    Sys.S = [Sys.S,S];
end
n = prod(2*Sys.S+1);
% each of the spins has a non zero B22 and B20 
B2 = [B22 0 B20 0 0]; % B(k=2,q) with q = +2,+1,0,-1,-2
 
Sys.B2 = [B2];
for i = 2:N_electrons
    Sys.B2 = [Sys.B2;B2];
end
J = 1.46*meV;
ee = [];ee1=[];
for i = 2:N_electrons
    ee1 = [1,zeros(1,i-2),ee1];
    %ee = [eye(3);zeros(3,3*(i-2));ee];
end
%Sys.ee = J.*ee;
if N_electrons>1;  Sys.J = ee1.*J; end

H=sham(Sys,[0,0,0],'sparse');
[Vecs,EE] = eig(full(H),'vector');





%-------------------------------------------------------------------------%
Sys1.S=Sys.S;
A0 = sparse(n,n);
A{1} = speye(n,n);
%Form Sparse A
A{2} = stev(Sys.S,2,0,1,'sparse');
A{3} = stev(Sys.S,2,2,1,'sparse');
for i = 2:N_electrons       
    A{2} = A{2} + stev(Sys.S,2,0,i,'sparse');
    A{3} = A{3} + stev(Sys.S,2,2,i,'sparse');
end   
l = length(A);
%Sys = rmfield(Sys,'J');
if N_electrons>1
    %Sys1.ee = ee; 
    Sys1.J = ee1;
    A{l+1} = eeint(Sys1,1:N_electrons,'sparse'); 
end

EE = EE(1:24);
assignin('base','EE',EE)
assignin('base','A',A)
assignin('base','A0',A0)
assignin('base',"H",H)
Sys.J = J.*ee1;
