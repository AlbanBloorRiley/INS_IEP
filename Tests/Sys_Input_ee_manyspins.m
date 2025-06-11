function ok = Sys_Input_ee_manyspins()
% Test for correct ordering of spin pairs

% Generate spin system with 3 coupled spins
N = 3;
Sys1.S = 1/2*ones(1,N);
pairs = nchoosek(N,2);
Sys1.J = rand(pairs,1);

% Add 3 more spins, without coupling to the first 3
Sys2.S = [Sys1.S,Sys1.S];
comb1 = nchoosek(1:N,2);
comb2 = nchoosek(1:(2*N),2);
pairs2 = nchoosek(2*N,2);
Sys2.J = zeros(pairs2,1);

% Copy interaction parameters
for n =1:pairs
  ind = logical((comb2(:,1) == comb1(n,1)) .* (comb2(:,2) == comb1(n,2)));
  Sys2.J(ind) = Sys1.J(n);
end

% Generat e-e interaction Hamiltonians
H1 = ham_ee(Sys1);
H2 = ham_ee(Sys2);

Vary1 = Sys1;
Vary2 = Sys2;

[A1,A01,x1] = Sys_Input(Sys1,Vary1);
[A2,A02,x2] = Sys_Input(Sys2,Vary2);


ok(1) = isempty(setdiff(A1(:),A2(:)));
ok(2) = isempty(setdiff(A01(:),A02(:)));
ok(3) = isempty(setdiff(x1(:),x2(:)));




% Determine whether there are elements that are different
% valdiff = setdiff(H1(:),H2(:));
% 
% ok = isempty(valdiff);