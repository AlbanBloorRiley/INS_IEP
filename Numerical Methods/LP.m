function [CurrentLoop,OutputLineLength] = LP(obj_fun,x0,~,params,RecordIterates)
NIter = 0;   x = x0;      OutputLineLength =[]; alpha = 0;
if RecordIterates, CurrentLoop.Iterates = x0; end
% lambda = params.method.Regularisation;
if isempty(params.method.ScalingMatrix); Scaling = false; else
    Scaling = true;
end

% Calculate residual, Jacobian and Hessian of R
[X.F,X.R,X.J] = obj_fun(x, params.method.constants); FuncCount = 1;
% Calculate deflation operators
% Print 0th iteration
if params.method.Verbose
    OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d; alpha = %d; \n', NIter,X.F,norm(X.J'*X.R),0);
end
[stop,CurrentLoop.ConvergenceFlag] = ismin(X.F, inf, X.J'*X.R, NIter, params.convergence);
% Main Loop
while stop == false


    % Print out convergence info
    if params.method.Verbose
        if NIter>0
            fprintf(repmat('\b',1,OutputLineLength))
        end
        OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d; alpha = %d; \n', NIter,X.F,norm(X.J'*X.R),alpha);
    end


    % A = FormA(x,params.method.constants.A,params.method.constants.A0);
    [QFull,DFull] = eig(full(FormA(x,params.method.constants.A,params.method.constants.A0)),'vector');
   
    if length(params.method.constants.ev)<length(DFull)
        C = (DFull'-params.method.constants.ev).^2;
        pairs = matchpairs(C,1e10);
        D = DFull(pairs(:,2));
        Q = QFull(:,pairs(:,2));
        params.method.constants.ev = params.method.constants.ev((pairs(:,1)));
    else
        D = DFull; Q = QFull;
    end


    Z = Q*diag(params.method.constants.ev)*Q';
    for i = 1:length(params.method.constants.A)
        b(i,1) = trace((Z-params.method.constants.A0)'*params.method.constants.A{i});
    end
    xprev = x;
    x = params.method.ScalingMatrix*b;
    
    NIter = NIter + 1;
    % Calculate residual, Jacobian of R
    [X.F,X.R,X.J] = obj_fun(x, params.method.constants); FuncCount = FuncCount +1;

        [stop, CurrentLoop.ConvergenceFlag] = ismin(X.F, x-xprev, X.J'*X.R, NIter, params.convergence);
    % Save iterates for plotting
    if RecordIterates; CurrentLoop.Iterates = [CurrentLoop.Iterates, x]; end
end
CurrentLoop.DeflatedPoint = x;
CurrentLoop.ErrorAtDeflatedPoint = obj_fun(x, params.method.constants);
CurrentLoop.NIter = NIter;
CurrentLoop.FuncCount = FuncCount;
end

