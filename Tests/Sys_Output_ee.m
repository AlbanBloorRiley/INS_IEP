function ok = Sys_Output_ee()

clear all
Sys.S = [7/2 7/2 7/2];   
Sys.B2 = rand(3,5);      
Sys.B4 =  rand(3,9);
Sys.ee = rand(3,3);

Vary.B2 = rand(3,5)<0.5;
Vary.B4 = rand(3,9)<0.5;
Vary.ee = rand(3,3)<0.5;

[~,~,x,Ops,SysFixed] = Sys_Input(Sys,Vary);
SysOut = Sys_Output(x,Ops,SysFixed);

ok(1) = isequal(Sys,SysOut);
ok(2) = areequal(ham(SysOut,[0,0,0]),ham(Sys,[0,0,0]),1e-10,'abs');