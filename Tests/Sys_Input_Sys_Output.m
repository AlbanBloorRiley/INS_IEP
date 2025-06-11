function ok = Sys_Input_Sys_Output()

Sys.S = [7/2 7/2];   
Sys.B2 = rand(2,5);      
Sys.B4 =  rand(2,9);
Sys.J = rand(1,1);

Vary.B2 = rand(2,5)<0.5;
Vary.B4 = rand(2,9)<0.5;
Vary.J = rand(1,1)<0.5;

[~,~,x,Ops,SysFixed] = Sys_Input(Sys,Vary);
SysOut = Sys_Output(x,Ops,SysFixed);

ok = isequal(Sys,SysOut);