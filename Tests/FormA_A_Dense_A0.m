function ok = FormA_A_Dense_A0()
% Test FormA when Ai are dense

N = 1e3;   n = 1e1;
x = rand(1,n);
A0 = rand(N);
Ax = A0;
for i = 1:n
    Ai{i} = rand(N);
    Ax = Ax + Ai{i}*x(i);
end

A = FormA(x,Ai,A0); 
ok(1) = areequal(A,Ax,1e-14,'rel');
ok(2) = ~issparse(A);