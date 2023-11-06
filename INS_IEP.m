function [SysOut, NIter, Flags, Iterates, FinalError] = INS_IEP(Sys0,Vary,Exp,varargin)
% INS_IEP Inelastic Neutron Scattering Inverse Eigenvalue Problem
% SysOut =  INS_IEP(Sys) produces an easy spin syle Sys structure
% containing the parameters found at a local minima of the error bewtween
% eigenvalues these simulate and the prescribed experimental eigenvalues input.
%
% The structure Sys, must contain the spin vector Sys.S and the
% experimental eigenvalues Sys.ev.
%
% The structure should also indicate which Stevens operators are to be fit
% (using Sys.B2,..,Sys.B6 etc. as done in easyspin)
% and, if using a multiple spin system, which exchange terms to use -
% either usung Sys.J xor Sys.ee. To indicate which are to be used make sure
% that the appropriate entry in the structure is non-zero, if using an initial guess please
% use this value. Make sure that no entries in the same field are
% identically the same value, unless you wish them to be pinned - in which
% case this is done by making the desired entries identically the same
% value. All values must be given in MH.
%
%[SysOut, NIter, Flags, Iters, FinalError] = INS_IEP(Sys) also produces
%cell structures containing the Number of Iterations till 'convergence';
% the type of 'convergence' achieved; all the iterations made until
% 'convergence'; the finall error at 'convergence'.
%
%Additional options: please pass these parameters as named optional
%arguments
%
%They are listed below as a tuple of
%"Name of Parameter" - Defalut Value - description:
%
%'NMinima' - (1) - The number of local minima you wish to find, the output
%Sys structures are given in an array, other outputs are in cell structures
%'Method' - {('Newton'),'GaussNewtonT1', 'GaussNewtonT1'} - The optimisation method desired
%'tol' - (1e-3) - The tolerence required for convergence
%'MaxIter' - (10000) - The integer value  of Maximum iterations per deflation,
%'Linesearch' - {('No','Basic', 'Armijo, 'Quadratic'} - line search method,
%'p' - (2) - The deflation exponent, an integer value or 'exp',
%'Eigensolver' - {'eig','eigs'} - The eigensolver desired, defalut
%chosed addaptively
%'A0' - zeros(n) - A single matrix shift to the hamiltonian, for example
%the Hamiltonian given by previously calculated Stevens operators
%'UseInitialGuess' - {true, (false)} - Set to true to use the input initial
%guess' of parameters,
%'c1' - 1e-4 - the armijo  line search parameter
%'alphaint' - 1 - Initial value of the line search parameter each iteration
%'tau' - 0.5 - The value of the decrease in the line search parameter.
%'minalpha' - 1e-4 - The minimu value of the line search parameter.
%deflatelinesearch' - true - Use deflated linesearch conditions or not.
%'theta' - 1 - The value of the shift in the deflation operator.
%'SysFound' - {([]),SysOut} - An array of structures of previously found local minima.
if nargin <3
    error('Sys0, Vary and Exp are required inputs')
end
if ~isfield(Exp,"ev")
    error("Please provide the experimental eigenvalues, using Exp.ev")
end

[A,A0,scale_x,Ops,SysFixed] = Sys_Input(Sys0,Vary);

% x0 = scale_x;
%  [A1,dint,ev1,Ops1] = SysInput(Sys0);
%     defaultDeflate = "F";




defaultMethod ='Newton';
defaultNMinima = 1;
if length(A{1})<1000
    defaultEigensolver = 'eig';
else
    defaultEigensolver = 'eigs';
end
defaultVerbose = false;
defaultScaled = true;


defaultLinesearch = 'Armijo';
defaultDeflateLinesearch = true;
defaultC1 = 1e-4;
defaultTau = 0.5;
defaultAlpha0 = 1;
defaultMinalpha = 1e-4;

defaultTheta = 2;
defaultSigma = 1;
defaultSingleShift = false;
defaultEpsilon = 0.1;


defaultGradientTolerance = 1e-2;
defaultStepTolerance = 1e-8;
defaultObjectiveTolerance = 1e-3;
defaultMaxIter = 10000;

defaultIEPType = 'Difference';
defaultSysFound = [];


% defaultA0 = sparse(length(A{1}),length(A{1}));

IP = inputParser;
addRequired(IP,'Sys0_In')
addRequired(IP,'Vary')
addRequired(IP,'Exp')

addParameter(IP,'Method',defaultMethod)
addParameter(IP,'NMinima',defaultNMinima)
addParameter(IP,'Verbose',defaultVerbose)
addParameter(IP,'Scaled',defaultScaled)

addParameter(IP,'Linesearch',defaultLinesearch)
addParameter(IP,'deflatelinesearch',defaultDeflateLinesearch)
addParameter(IP,'c1',defaultC1)
addParameter(IP,'tau',defaultTau)
addParameter(IP,'alpha0',defaultAlpha0)
addParameter(IP,'minalpha',defaultMinalpha)


addParameter(IP,'theta',defaultTheta)
addParameter(IP,'sigma',defaultSigma)
addParameter(IP,'SingleShift',defaultSingleShift)
addParameter(IP,'epsilon',defaultEpsilon)


addParameter(IP,'ObjectiveTolerance',defaultObjectiveTolerance)
addParameter(IP,'GradientTolerance',defaultGradientTolerance)
addParameter(IP,'StepTolerance',defaultStepTolerance)
addParameter(IP,'MaxIter',defaultMaxIter)

addParameter(IP,'IEPType',defaultIEPType)
addParameter(IP,'SysFound',defaultSysFound)
addParameter(IP,'Eigensolver',defaultEigensolver)


IP.parse(Sys0,Vary,Exp,varargin{:})
res = IP.Results;


%Method Parameters
params.method.StepMethod = res.Method;
params.method.MinimaRequested = res.NMinima;
params.method.Verbose = res.Verbose;
params.method.Scaled = res.Scaled;

%Line search Parameters
params.linesearch.method = res.Linesearch;
params.linesearch.deflate = res.deflatelinesearch;
params.linesearch.c1 = res.c1;
params.linesearch.tau = res.tau;
params.linesearch.alpha0 = res.alpha0;
params.linesearch.minalpha = res.minalpha;

%Deflation Parameters
params.deflation.theta = res.theta;
params.deflation.sigma = res.sigma;
params.deflation.singleshift = res.SingleShift;
params.deflation.epsilon = res.epsilon;


%Convergence Parameters
params.convergence.ObjectiveTolerance = res.ObjectiveTolerance;
params.convergence.GradientTolerance = res.GradientTolerance;
params.convergence.StepTolerance = res.StepTolerance;
params.convergence.MaximumIterations = res.MaxIter;

% %INS IEP type
params.IEP.Eigensolver = res.Eigensolver;

% params.IEP.type = res.IEPType;


if res.Scaled
    for i =1:length(scale_x)
        A{i}=scale_x(i)*A{i};
    end
    x0=ones(length(scale_x),1);
else
    x0=scale_x;
end



% if ~res.UseInitialGuess
%     x0 = ones(length(A),1);
% end

Y=[];
for i= 1:length(res.SysFound)
    [~,~,Y(:,i)] = Sys_Input(res.SysFound(i),Vary);
end
if size(Exp.ev,1)<size(Exp.ev,2)
    constants.ev = Exp.ev';
else
    constants.ev = Exp.ev;
end


if res.IEPType == "Classic"
    obj_fun = @IEP_Evaluate_full;
    A{end+1} = speye(size(A{1}));
    x0(end+1)=1;
    if isfield(Exp,'evsd')
        if ~all(Exp.evsd)
            warning('If standard deviation is given for experimental eigenvalues, any zero values are set to 1.')
            Exp.evsd(~Exp.evsd) = 1;
        end
        constants.eigenvalueSD = Exp.evsd;
    else
        if isfield(Exp,'evdsd')
            warning('If IEP type is set to ''Classic'' then Exp.evdsd is not used. Please use Exp.evsd.')
        end
        constants.eigenvalueSD = ones(length(Exp.ev),1);
    end

elseif res.IEPType == "Difference"
    obj_fun = @IEP_Evaluate_diff;

    if isfield(Exp,'evdsd')
        if ~all(Exp.evdsd)
            warning('If standard deviation is given for experimental eigenvalue differences, any zero values are set to 1.')
            Exp.evsd(~Exp.evsd) = 1;
        end                     
        constants.eigenvalueDifferenceSD = Exp.evdsd;
    else
        if isfield(Exp,'evsd')
            warning('If IEP type is set to ''Difference'' then Exp.evsd is not used. Please use Exp.evdsd.')
        end
        constants.eigenvalueDifferenceSD = ones(length(Exp.ev)-1,1);
    end

else
    error('Please use either the ''Classic'' or ''Difference'' options for IEPType')
end
constants.A = A;
constants.A0 = A0;
constants.ED = res.Eigensolver;
params.method.constants = constants;



% Set line search objective/merit function and derivatives - phi/gradphi
switch res.Method
    case "Newton"
        if res.deflatelinesearch
%             params.linesearch.phi_0=@(obj_at_x,Mu,Rx,Jx) Mu^2*(Jx'*Rx)'*(Jx'*Rx);
            params.linesearch.phi = @(objective_function,x,constants,PreviousMinima,DeflationParameters) Deflated_Gradient(objective_function,x,constants,PreviousMinima,DeflationParameters);
%             params.linesearch.gradphi = @(Rx,Jx,Mu,gradMu,S) (Jx'*Rx*gradMu+Mu*(Jx'*Jx+S))*(2*Mu*Jx'*Rx);
            params.linesearch.gradphi = @(Rx,Jx,Mu,gradMu,S) 2*Mu*gradMu'*(Jx'*Rx)'*(Jx'*Rx)+ 2*Mu^2*(Jx'*Jx+S)*(Jx'*Rx);

        else
%             params.linesearch.phi_0=@(obj_at_x,Mu,Rx,Jx)obj_at_x;
            params.linesearch.phi = @(objective_function,x,constants,PreviousMinima,DeflationParameters)  objective_function(x,constants);
            params.linesearch.gradphi = @(Rx,Jx,Mu,gradMu,S) 2*(Jx'*Rx);

        end

        Fun= str2func('Deflated_Newton_Method');

    case "Gauss-NewtonT1"
        if res.deflatelinesearch
            if isfield(varargin{1},"deflatelinesearch")
                warning("There is no deflated line search method for the Type 1 Gauss-Newton method. Using an undeflated line search.")
            end
        end
%         params.linesearch.phi_0=@(obj_at_x,Mu,Rx,Jx)obj_at_x;
        params.linesearch.phi = @(objective_function,x,constants,PreviousMinima,DeflationParameters) objective_function(x,constants);
        params.linesearch.gradphi = @(Rx,Jx,Mu,gradMu,S) 2*(Jx'*Rx);

        Fun= str2func('Deflated_GaussNewton_Type1_Method');

    case "Gauss-NewtonT2"
        if res.deflatelinesearch
%             params.linesearch.phi_0=@(obj_at_x,Mu,Rx,Jx)Mu^2*obj_at_x;
            params.linesearch.phi = @(objective_function,x,constants,PreviousMinima,DeflationParameters) deflation(PreviousMinima,x,DeflationParameters)^2*objective_function(x,constants);
            params.linesearch.gradphi = @(Rx,Jx,Mu,gradMu,S) 2*(Jx'*Rx);
        else
%             params.linesearch.phi_0=@(obj_at_x,Mu,Rx,Jx)obj_at_x;
            params.linesearch.phi = @(objective_function,x,constants,PreviousMinima,DeflationParameters) objective_function(x,constants);
            params.linesearch.gradphi = @(Rx,Jx,Mu,gradMu,S) 2*(Jx'*Rx);
        end

        Fun= str2func('Deflated_GaussNewton_Type2_Method');

    otherwise
        error('Please select a valid method - Newton/Gauss-NewtonT1/Gauss-NewtonT2')
end
%[x,NIter,flag,obj_at_x,nbytes,xs]
% [SysOut, NIter, Flags, Iterates, FinalError] 
j=0;tic;
for i=length(res.SysFound)+1:res.NMinima
    [Y(:,i),NIter(i),Flags{i},FinalError{i},nbytes,Iterates{i}] = Fun(obj_fun,x0,Y,params);
    if Flags{i} ~="Max Iterations reached"
        j = j+1;
        disp(['Currently found ',num2str(j),' minima in ',num2str(i),' iterations.'])
    else
        disp(['Currently found ',num2str(j),' minima in ',num2str(i),' iterations.'])
    end
end; toc


if res.IEPType == "Classic"
    Y=Y(1:end-1,:);
end
if res.Scaled
    SysOut = Sys_Output(Y.*scale_x,Ops,SysFixed);
else
    SysOut = Sys_Output(Y,Ops,SysFixed);
end

end

function f_out = Deflated_Gradient(objective_function,x,constants,PreviousMinima,DeflationParameters)
[~,R,J] = objective_function(x,constants); 
f_out = deflation(PreviousMinima,x,DeflationParameters)^2*(J'*R)'*(J'*R);
end