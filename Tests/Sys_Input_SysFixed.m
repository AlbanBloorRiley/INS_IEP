function ok = Sys_Input_SysFixed()

Sys.S = [7/2 7/2];   
Sys.B2 = rand(2,5);      
Sys.B4 =  rand(2,9);
Sys.J = rand(1,1);

Vary.B2 = rand(2,5)<0.5;
Vary.B4 = rand(2,9)<0.5;
Vary.J = rand(1,1)<0.5;

[~,~,~,~,SysFixed] = Sys_Input(Sys,Vary);

Sys.B2(Vary.B2) = 0;
Sys.B4(Vary.B4) = 0;
Sys.J(Vary.J) = 0;

ok(1) = isequal(SysFixed,Sys);
ok(2) = areequal(ham(SysFixed,[0,0,0]),ham(Sys,[0,0,0]),1e-10,'abs');
