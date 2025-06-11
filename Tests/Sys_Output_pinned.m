function ok = Sys_Output_pinned()

Sys.S = [7/2 7/2];   
Sys.B2 = rand(1,5).*[1;1];      
Sys.B4 =  rand(1,9).*[1;1];
Sys.J = rand(1,1);

Vary.B2 = rand(1,5).*[1;1]<0.5;
Vary.B4 = rand(1,9).*[1;1]<0.5;
Vary.J = rand(1,1)<0.5;

[~,~,x,Ops,SysFixed] = Sys_Input(Sys,Vary);
SysOut = Sys_Output(x,Ops,SysFixed);

ok(1) = true;
for i = 1:length(Ops.B2)
    if Ops.B2{i}(1,:)~= Ops.B2{i}(2,:)
        ok(1) = false;
        break
    end
end
ok(2) = true;
for i = 1:length(Ops.B4)
    if Ops.B4{i}(1,:)~= Ops.B4{i}(2,:)
        ok(2) = false;
        break
    end
end

ok(3) = isequal(Sys,SysOut);