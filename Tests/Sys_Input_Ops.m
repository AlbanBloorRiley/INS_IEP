function ok = Sys_Input_Ops()

Sys.S = [7/2 7/2];   
Sys.B2 = rand(2,5);      
Sys.B4 =  rand(2,9);
Sys.ee = rand(1,3);

Vary1.B2 = rand(2,5)<0.5;
Vary1.B4 = rand(2,9)<0.5;
Vary1.ee = rand(1,3)<0.5;

[~,~,~,Ops] = Sys_Input(Sys,Vary1);

Vary2.B2 = zeros(size(Vary1.B2));
for i = 1:length(Ops.B2)
    Vary2.B2 = Vary2.B2 +Ops.B2{i};
end
ok(1) = areequal(double(Vary1.B2),Vary2.B2);

Vary2.B4 = zeros(size(Vary1.B4));
for i = 1:length(Ops.B4)
    Vary2.B4 = Vary2.B4 +Ops.B4{i};
end
ok(2) = areequal(double(Vary1.B4),Vary2.B4);

Vary2.ee = zeros(size(Vary1.ee));
for i = 1:length(Ops.ee)
    Vary2.ee = Vary2.ee +Ops.ee{i};
end
ok(3) = areequal(double(Vary1.ee),Vary2.ee);