function [CurrentLoop] = RGD_LP(obj_fun,x,~,params,RecordIterates)
NIter = 0;   constants = params.method.constants;
CurrentLoop.Iterates = x;
Fprev = inf; 
Binv = params.method.ScalingMatrix*FormBinv(constants.A);

% Calculate residual, Jacobian and Hessian of R
[X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = 1;
% CurrentLoop.Error = X.F;
if params.method.Verbose
    OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d\n', NIter,X.F,norm(X.J'*X.R));
    fprintf(repmat(' ',1,OutputLineLength))
end
[stop,CurrentLoop.ConvergenceFlag] = ismin(X.F,x ,inf,X.J'*X.F, NIter, params.convergence);
% Main Loop )
while stop == false
    p = - Binv*X.J'*X.R;
    xprev = x;

    x = xprev + p;
    NIter = NIter + 1;
    % Calculate residual, Jacobian of R
    [X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = FuncCount +1;


    if Fprev<X.F
        X.F;
    end
    Fprev = X.F;
    % [constants,x,Binv] = rescale(constants,x);
    if params.method.Verbose
        fprintf(repmat('\b',1,OutputLineLength))
        OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d \n', NIter,X.F,norm(X.J'*X.R));
    end
    [stop, CurrentLoop.ConvergenceFlag] = ismin(X.F,x, p, X.J'*X.F,NIter, params.convergence);
    pprev=p;
    % Save iterates for plotting
    if RecordIterates
        CurrentLoop.Iterates = [CurrentLoop.Iterates, x];
    end
end
CurrentLoop.DeflatedPoint = x;
CurrentLoop.ErrorAtDeflatedPoint = obj_fun(x, constants);
CurrentLoop.NIter = NIter;
CurrentLoop.FuncCount = FuncCount;
end


