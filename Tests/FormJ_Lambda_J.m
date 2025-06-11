function ok = FormJ_Lambda_J()
N = 10;
l = 4;
for i = 1:l
    X = rand(N)*1e2;
    A{i} = X'*X;
end
x0 = rand(l,1);
ev = rand(N,1);

[Q,~] = eig(FormA(x0,A));
J = FormJ_Lambda(Q,A);

fun = @(x)(eig(FormA(x,A)) - ev);
jac = jacobianest(fun,x0);
ok = areequal(J,jac,1e-8,'rel');
