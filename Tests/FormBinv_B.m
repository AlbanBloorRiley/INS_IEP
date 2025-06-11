function ok = FormBinv_B()

A{1} = [0 1; 1 2];
A{2} = [1 0; 0 1];
A{3} = [0 0; 1 2];
[B1inv, B1] = FormBinv(A);

B2 = [6     2     5;
      2     2     2;
      5     2     4];

B2inv = [ -2           -1            3;
          -1          0.5            1;
           3            1           -4];

ok(1) = areequal(B1,B2);
ok(2) = areequal(B1inv,B2inv, 1e-14, 'abs');
ok(3) = areequal(B1*B1inv, eye(3), 1e-14, 'abs');