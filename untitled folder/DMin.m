function [Y, NIter, Flags, Iters, FinalError] = DMin(Evaluate,xint,varargin)
%DMin Finds multiple minima of the problem using deflation.
%
% Y= DMin(Evaluate,xint,varargin)   Outputs an array of vectors that
% minimise the problem 'Evaluate', initialised as xint.
%    * Evaluate(x,constants) is a function handle that returns the current 
%    error of the least squares formulation, given the current guess, x, 
%    and optional parameters, constants. It must also optionally output the
%    Residual and the  Jacbian and Hessian of said residual of the least
%    squares problem: [F,R,J,H] = Evaluate(x,constants)
%
%[Y, NIter, Flags, Iters, FinalError] = DMin(Evaluate,xint,varargin) also 
% produces cell structures containing the Number of Iterations till 
% convergence; the type of 'convergence' achieved; all the iterations 
% made until convergence; the finall error at convergence.
%
%Additional options: please pass these parameters as named optional
%arguments
%
%They are listed below as a tuple of
%"Name of Parameter" - Defalut Value - description:
%
%'NMinima' - (1) - The number of local minima you wish to find, the output
%'Method' - {('Newton'),'GaussNewtonT1', 'GaussNewtonT1'} - The optimisation method desired
%'tol' - (1e-3) - The tolerence required for convergence
%'MaxIter' - (10000) - The integer value  of Maximum iterations per deflation,
%'Linesearch' - {('No','Basic', 'Armijo, 'Quadratic'} - line search method,
%'p' - (2) - The deflation exponent, an integer value or 'exp',
%'UseInitialGuess' - {true, (false)} - Set to true to use the input initial
%guess' of parameters,
%'UseInitialGuess' - {true, (false)} - Set to true to use the input initial
%guess' of parameters,
%'c1' - 1e-4 - the armijo  line search parameter
%'alphaint' - 1 - Initial value of the line search parameter each iteration
%'tau' - 0.5 - The value of the decrease in the line search parameter.
%'minalpha' - 1e-4 - The minimu value of the line search parameter.
%deflatelinesearch' - true - Use deflated linesearch conditions or not.
%'theta' - 1 - The value of the shift in the deflation operator.

defaultUseInitialGuess = false;
defaultNMinima = 1;
defaultLinesearch = 'No';
defaultTolerance = 1e-3;
defaultMaxIter = 10000;
defaultP = 2;
defaultC1 = 1e-4;
defaultAlphaInt = 1;
defaultTau = 0.5;
defaultMinalpha = 1e-2;
defaultTheta = 1;
defaultMethod = 'Gauss-NewtonT1';
defaultConstants = [];
defaultDeflateLinesearch = true;



IP = inputParser;
addRequired(IP,'Evaluate')
addRequired(IP,'xint')
addParameter(IP,'Method',defaultMethod)
addParameter(IP,'constants',defaultConstants)
addParameter(IP,'UseInitialGuess',defaultUseInitialGuess)
addParameter(IP,'NMinima',defaultNMinima)
addParameter(IP,'Linesearch',defaultLinesearch)
addParameter(IP,'Tolerance',defaultTolerance)
addParameter(IP,'MaxIter',defaultMaxIter)
addParameter(IP,'p',defaultP)
addParameter(IP,'c1',defaultC1)
addParameter(IP,'alphaint',defaultAlphaInt)
addParameter(IP,'tau',defaultTau)
addParameter(IP,'minalpha',defaultMinalpha)
addParameter(IP,'theta',defaultTheta)
addParameter(IP,'deflatelinesearch',defaultDeflateLinesearch)
IP.parse(Evaluate,xint,varargin{:})

res = IP.Results;

Y=[];

switch res.Method
    case "Newton"
        if res.deflatelinesearch
            f_0="M^2*(Jx'*Rx)'*(Jx'*Rx)";
            f = @(x,constants,Evaluate,Y,p)F_JDeflated(x,constants,Evaluate,Y,p);
            J_f = "(Jx'*Rx*gradM+M*(Jx'*Jx+S))*(2*M*Jx'*Rx)";
        else
            f_0="(Jx'*Rx)'*(Jx'*Rx)";
            f = @(x,constants,Evaluate,Y,p)Evaluate(x,constants);
            J_f = "(Jx'*Jx+S)*(Jx'*Rx)";
        end
        Fun= str2func('Deflated_Newton_Method');
    case "Gauss-NewtonT1"
        if res.deflatelinesearch
            f_0="M^2*(Jx'*Rx)'*(Jx'*Rx)";
            f = @(x,constants,Evaluate,Y,p)F_JDeflated(x,constants,Evaluate,Y,p);
            J_f = "2*(Jx'*Rx*gradM+M*(Jx'*Jx))*(2*M*Jx'*Rx)";
        else
            f_0="(Jx'*Rx)'*(Jx'*Rx)";
            f = @(x,constants,Evaluate,Y,p)Evaluate(x,constants);
            J_f = "2*(Jx'*Jx)*(Jx'*Rx)";
        end
        Fun= str2func('Deflated_GaussNewton_Type1_Method');
    case "Gauss-NewtonT2"
        if res.deflatelinesearch
            f_0="M^2*Fprev";
            f = @(x,constants,Evaluate,Y,p)F_RDeflated(x,constants,Evaluate,Y,p);
            J_f = "2*M*gradM + 2*M^2*Jx'*Rx";
        else
            f_0="Fprev";
            f = @(x,constants,Evaluate,Y,p)Evaluate(x,constants);
            J_f = "2*Jx'*Rx";
        end
        Fun= str2func('Deflated_GaussNewton_Type1_Method');
    otherwise
        error('Please select a valid method - Newton/Gauss-NewtonT1/Gauss-NewtonT2')
end

for  i=1:res.NMinima
    [Y(:,i),NIter(i),Flags{i},Iters{i},FinalError{i}] = ...
        Fun(Evaluate,xint,res.constants,res.Tolerance,...
        res.MaxIter,Y,res.p,res.Linesearch,res.c1,res.alphaint,...
        res.tau,res.minalpha,res.theta,f,f_0,J_f);
    disp(['Currently found ',num2str(i),' minima'])
end
end



function f_out = F_JDeflated(d,constants,Evaluate,Y,p)
M = deflation(Y,d,p);
[~,R,J] = Evaluate(d,constants); f_out = M^2*(J'*R)'*(J'*R);
end
function f_out = F_RDeflated(d,constants,Evaluate,Y,p)
M = deflation(Y,d,p);
f_out = M^2*Evaluate(d,constants);
end