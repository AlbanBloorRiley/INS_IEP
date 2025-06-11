function ok = Sys_Output_multiplesolutions()

N = 15;
Vary.B2 = rand(1,5).*[1;1]<0.5;
Vary.B4 = rand(1,9).*[1;1]<0.5;
Vary.J = rand(1,1)<0.5;


for i = 1:N
Sys(i).S = [7/2 7/2];   
Sys(i).B2 = rand(1,5).*Vary.B2;      
Sys(i).B4 =  rand(1,9).*Vary.B4;
Sys(i).J = rand(1,1).*Vary.J;
end

[~,~,x,Ops,SysFixed] = Sys_Input(Sys(1),Vary);

X=x;
for i = 2:N
    [~,~,x] = Sys_Input(Sys(i),Vary);
    X = [X,x];
end
SysOut = Sys_Output(X,Ops,SysFixed);
ok = true;
for i = 1:N
    if ~ areequal(ham(SysOut(i),[0,0,0]),ham(Sys(i),[0,0,0]),1e-10,'abs')
        ok = false;
    end
end
