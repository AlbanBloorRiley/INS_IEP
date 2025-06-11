function ok = SysInput_ee_full()
% Tests Sys_input computation of basis matrices, for electron-electron 
% interaction terms.

ee1 = [1 1 1];
ee2 = [2 2 2];
ee3 = [5 5 5];

Sys1.S = [1/2 1/2 1/2];
Sys1.ee = [ee1; ee2; ee3];
Sys2.S = Sys1.S;
Sys2.ee = [diag(ee1);diag(ee2);diag(ee3)];
Vary1 = Sys1;
Vary2 = Sys2;

[A1,A01,x1] = Sys_Input(Sys1,Vary1);
[A2,A02,x2] = Sys_Input(Sys2,Vary2);

ok(1) = areequal(cell2mat(A1),cell2mat(A2),1e-10,'abs');
ok(2) = areequal(A01,A02,1e-10,'abs');
ok(3) = areequal(x1,x2,1e-10,'abs');


