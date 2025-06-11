function ok = FormA_A_Sparse()
% Test FormA when Ai are sparse

N = 1e3;   n = 1e1;
x = rand(1,n);
Ax = zeros(N);

for i = 1:n
    Ai{i} = sprand(N,N,0.1);
    Ax = Ax + Ai{i}*x(i);
end
 A = FormA(x,Ai);
ok(1) = areequal(A,Ax,1e-14,'rel');
ok(2) = issparse(A);