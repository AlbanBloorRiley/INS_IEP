function ok = FormH_Lambda_H()

N = 10;
l = 4;
for i = 1:l
    X = rand(N);
    A{i} = X'*X;
end
x0 = rand(l,1);
ev = rand(N,1);

[Q,D] = eig(FormA(x0,A));
H = FormH_Lambda(Q,A,D,ev);


for i = 1:N
    fun = @(x)(itheigenvalue(FormA(x,A),i)  - ev(i));
    hes(:,:,i) = hessian(fun,x0);
end

ok = areequal(H,hes,1e-7,'rel');

end


function e = itheigenvalue(A,i)
    E = eig(A);
    e = E(i);
end