function ok = Sys_Input_twoelectrons()

Sys.S = [7/2 7/2];   
Sys.B2 = rand(2,5);      
Sys.B4 =  rand(2,9);
Sys.J = 0;

Vary.B2 = rand(2,5)<0.5;
Vary.B4 = rand(2,9)<0.5;
Vary.J = 0;

[A,A0,x] = Sys_Input(Sys,Vary);

ok = areequal(FormA(x,A,A0),ham(Sys,[0,0,0]),1e-10,'abs');
