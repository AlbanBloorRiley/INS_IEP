function ok = Linesearch_Armijo()
l = 5; 

A = rand(l);
B = rand(l);
C = rand(l,1);
f = @(x,constants) norm(A*x + x'*B*x + C);

x = rand(l,1)*1e1;



params.method.constants = [];
params.deflation = [];

linesearch.phi = @(objective_function,x,constants,~,~) objective_function(x,constants);
linesearch.gradphi = @(X) 2*(X.J'*X.R);
linesearch.alpha0 = 1;
linesearch.method = "Armijo";
linesearch.c1 = 1e-4;
linesearch.tau = 0.5;
linesearch.minalpha = 1e-18;


X.F = f(x,[]);
X.J = A + 2*B*x;
X.R = A*x + x'*B*x + C;

%gradient descent 
p = -X.J'*X.R;

alpha = Linesearch(x,p,[],params,linesearch,f,X);

ok =  f(x+alpha*p) < f(x);
