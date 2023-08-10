
Sys = Mn12_Spin_Sys_3(2,2);
[SysOut, NIter, Flags, Iters, FinalError] =INS_IEP(Sys,'method','Gauss-NewtonT1','UseInitialGuess',true,'NMinima',4,'Linesearch','No','minalpha',1e-2);
%%
constants.ev = Exp.ev;

[constants.A,constants.A0,scale_x,Ops] = Sys_Input(Sys,Vary);
% constants.A0 = sparse(length(constants.A{1}),length(constants.A{1}));
constants.ED ='eig';

if Opt.Scaled
    for i =1:length(scale_x)
        constants.A{i}=scale_x(i)*constants.A{i};
    end
end

obj_fun = @(x)IEP_Evaluate_diff(x,constants);


