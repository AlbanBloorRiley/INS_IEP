function [test,flag] = ismin(f,x,p,gradf,NIter,ConvergenceParams,varargin)
%ISMIN Tests stopping criteria at current iteration
%
% stop = ismin(f,p,x,gradf,NIter,ConvergenceParams) outputs true if
% stopping critera reached, and false if else. f is the current value of
% the function evaluation, p is the current update step, x is the current
% iteration, gradf is the gradient of f at the current iteration, NIter is
% the iteration number, and ConvergenceParams is a structure containing the
% tolerences for each stopping criterion.
%
% [stop,flag] = ismin(f,p,gradf,NIter,ConvergenceParams) also outputs
% which stopping criteria was reached.

test = false;
flag = "";
if norm(p)<ConvergenceParams.StepTolerance
    flag = 'Step Size below tolerance';
    test=true;
elseif isnan(f)||isinf(f)
    test = true ;
    flag = 'NaN/Inf';
elseif f <ConvergenceParams.FunctionTolerance
    flag = 'Objective less than tolerance';
    test = true ;
elseif norm(gradf)<ConvergenceParams.GradientTolerance
    flag = 'Gradient less than tolerance';
    test=true;
elseif NIter>=ConvergenceParams.MaximumIterations
    flag ='Max Iterations reached';
    test = true  ;
elseif nargin> 6 && (any(isnan(varargin{1}))||any(isinf(varargin{1})))
    flag = "Deflation operator is Nan/Inf";
    test = true;
    elseif isfield(ConvergenceParams,'RelativeStepTolerance')&& (norm(p)/norm(x))<ConvergenceParams.RelativeStepTolerance
    flag = 'Relative Step Size below tolerance';
    test=true;
elseif isfield(ConvergenceParams,'StepDifferenceTolerance')&& abs(norm(p)-norm(varargin{2}))<ConvergenceParams.StepDifferenceTolerance
    flag = 'Step Size Difference below tolerance';
    test=true;
end
end
