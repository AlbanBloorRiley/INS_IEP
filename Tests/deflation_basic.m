function ok = deflation_basic()

x = rand(10,1)*10;
Y = rand(10,100)*10;
warning('off')
ok(1) = isinf(deflation(x,x));
warning('off')
ok(2) = ~isinf(deflation(x,x+1));
ok(3) = deflation(Y(:,1),x)>1;
ok(4) = deflation(Y(:,1),x+1e5)-1<1e-4;
ok(5) = deflation(Y(:,1),x) < deflation(Y,x);


