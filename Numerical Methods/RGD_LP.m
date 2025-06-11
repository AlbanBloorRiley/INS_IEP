function [CurrentLoop] = RGD_LP(obj_fun,x,DeflatedPts,params,RecordIterates)
NIter = 0;   constants = params.method.constants;
CurrentLoop.Iterates = x;
Fprev = inf; 
Binv = params.method.ScalingMatrix*FormBinv(constants.A);

% Calculate residual, Jacobian and Hessian of R
[X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = 1;
% Calculate updated deflation operators
[X.Mu,X.gradMu] = deflation(DeflatedPts, x, params.deflation);
% CurrentLoop.Error = X.F;
if params.method.Verbose
    OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d\n', NIter,X.F,norm(X.J'*X.R));
    fprintf(repmat(' ',1,OutputLineLength))
end
[stop,CurrentLoop.ConvergenceFlag] = ismin(X.F,x ,inf,X.J'*X.F, NIter, params.convergence);
% Main Loop )
while stop == false
    p = - params.method.ScalingMatrix*Binv*X.J'*X.R;

     pTgradNu = dot(X.gradMu,p)/X.Mu;

    alpha = params.linesearch.merit.alpha0;
    if pTgradNu > 1 % Zone 3
        % Could use a linesearch on mu here. Not necessary in practice.
        alpha = alpha/(1-pTgradNu);
    elseif pTgradNu >= params.deflation.epsilon % Zone 2
        % Perform a deflated step (no line search, but possibly alpha0 damping)
        alpha = alpha/(1-pTgradNu);
    end
    % xprev = x;
    x = x + alpha*p;
    NIter = NIter + 1;
    % Calculate residual, Jacobian of R
    [X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = FuncCount +1;
    % Calculate updated deflation operators
    [X.Mu,X.gradMu] = deflation(DeflatedPts, x, params.deflation);

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
    % pprev=p;
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


