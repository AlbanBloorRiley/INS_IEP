function ok = Sys_Input_ZF()
% Tests Sys_input computation of basis matrices, for zero field splitting,
% with some fixed values.

Sys.S = 7/2;   
Sys.B2 = (rand(1,5)-0.5)*1e5;      
Sys.B4 =  (rand(1,9)-0.5)*1e5;

Vary.B2 = rand(1,5)<0.5;
Vary.B4 = rand(1,9)<0.5;


[A,A0,x] = Sys_Input(Sys,Vary);

% H = ham(Sys,[0,0,0]);
% tol = 1e-1;
% all(all(FormA(x,A,A0) - H < tol))

ok = isequal(FormA(x,A,A0),ham(Sys,[0,0,0]));
