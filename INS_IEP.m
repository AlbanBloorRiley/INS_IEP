function [SysOut,Opt,params,obj_fun,Iterations] = INS_IEP(Sys0,varySys,Exp,varargin)
% INS_IEP Inelastic Neutron Scattering Inverse Eigenvalue Problem
% SysOut =  INS_IEP(Sys0,Vary,Exp) produces an easyspin style Sys structure
% containing the parameters that minimise difference in between the
% eigenvalues simulated and the prescribed experimental eigenvalues input.
% The function models closely the style and use of easyspin's esfit.
%
% The required inputs are Sys, Vary and Exp. Sys is used to specify the
% model being used to simulate the system, using Spin, Zero Field Splittings and
% Electron-electron spin-spin couplings, and gives the initial values for
% all parameters of the model. Vary specifies which paramaters are to be
% varied, and Exp contains the experimental data. All inputs conform to the
% easyspin syntax unless explicitly stated.
%
% Defining the Sys structure:
% Spin - As in easy spin is defined Sys.S = S (1xN array of real)
% Zero Field Splittings - These are given using the extended Stevens
% operators, using easyspin's syntax for "High-order zero-field splittings"
% Sys.B2,..,Sys.B12. Not all have to be used. Any entries with the same
% value within the same column (representing the same operator across
% multiple spin centres) will be pinned.
% Electron-electron spin-spin couplings -  given by either using Sys.J
% (1xN array of real) or Sys.ee (Nx3 or 3Nx3 array of real) as in easyspin.
% Note thot only one of these may be used at one time Note also that any
% entries of J or diagonals of ee that are equal will be pinned, also any
% symmetric ellements of ee that are equal to each others reciprical.
%
% Specifying Vary:
% Vary should contain all the same fields of Sys (not including S), any
% entry that is non-zero specifies that it is a parameter to be varied,
% using the initial value given in Sys. Note that parameters that are given
% by Sys but are not to be varied will be used as fixed values.
%
% Experimental data:
% Exp should contain the experimentally calculated eigenvalues saved as
% Exp.ev. It can optionally also give the associated uncertainty of the
% eigenvalues as Exp.evsd or uncertainty of the difference in eigenvalues
% as Exp.evdsd. These two options come from the two different formulations
% of the Inverse Eigenvalue Problem - where the residual is given as the
% differnce in eigenvalues, or the difference in the difference between
% adjacent eigenvalues. The specific formulation to be used is specified by
% the option 'IEPType' as either 'Classic' or 'Difference' respectivly.
%
%
% OPTIONAL PARAMETERS
% INS options:
% IEPType - Defines which formulation of the IEP to use - [ "Classic" |
%        {"Difference"} ]
% SysFound - array of Sys structures of previously found minimising systems
%  - [ {[]}, SysFound]
% GroundStateFound - Value of groundstate found if not already specified by
%        SysFound - [ [] | scalar ]
% Main method options:
% Method - The optimisation method to use - [ {"Good_GN"} | "Bad_GN" |
%            "Newton" | "RGD_LP" ]
% NDeflations - The number of local minima you wish to find -
%  [ {1} | positive integer ]
% MaxNonMinima - number of failed deflations before quiting method - [ {inf}
%         | positive integer ]
% Solver Options:
% Regularisation - The value regularising parameter - [ {0} | 
%          positive scalar ]
% LinearSolver - The linear solver used - [ {"mldivide"} | "lsqminnorm" ]
% Scaled - specifies if the parameters are scaled to the size of the initil
%        values - [ true | {false} ]
% Eigensolver - Specifies which eigensolver to use - [ "eig" | "eigs" ] 
%        (default is decided based on the size of the matrices used)
% Input/Output Options:
% Verbose - Output value of objective function and gradient at each
%             iteration - [ true | {false} ]
% RecordIterates - Option to not record all iterates if running for large
%        number of iterations - [ {true} | false ]
% SupressConvergedFlag - Stops outputing when method has converged -  
%        [ {false} | true ]
% EigsNotConvergedWarning - Specifies if a warning should be output of the
%        eigensolver fails to converge for some of the eigenvalues - 
%        [ {true} | false ]
% Lines earch options:
% Linesearch - line search method used - [ "No" | "Basic" | {"Armijo"} | 
%               "Quadratic" ]
% C1 - the armijo line search parameter - [ {1e-4} |
%          positive scalar in  [0,1] ]
% Tau - The value of the decrease in the line search parameter - [ {0.5} |
%          positive scalar in [0,1] ]
% alpha0 - Initial value of the line search parameter each iteration - [ {1}
%            | positive scalar ]
% Minalpha - The minimum value of the line search parameter - [ {1e-15 |
%          positive scalar in [0,1] ]
% DeflatedLinesearch - Line search on deflated objective for
%        Bad Gauss-Newton method - [ {false} | true ]
% Deflation options:
% Theta - The deflation exponent - [ {2} | positive integer | "exp" ]
% Sigma - The value of the shift of the deflation - [ {1} | positive scalar ]
% ShiftType - Deflation shift should only be applied once [ true |
%         {false} ]
% Epsilon - Tolerence on the application of deflation or line search - [
%         {0.01} | scalar ]
% Stopping criteria:
% StepTolerance - Stopping criterion based on difference of consecutive
%         steps - [ {1e-6} | positive scalar ]
% GradientTolerance - Stopping criterion based on gradient of function -
%         [ {0} | positive scalar ]
% FunctionTolerance - Stopping criterion based on value of objective function
%         function - [ {0} | positive scalar ]
% RelativeStepTolerance - Stopping criterion based on value of objective function
%         function - [ {0} | positive scalar ]
% MaxIter - The integer value  of Maximum iterations per deflation -
%            [ {1000} | positive integer ]
% ConvergenceFlags - List flags that are considered to mean convergence - [
%         {"Objective less than tolerance","Gradient less than tolerance",
%        "Step Size below tolerance","Merit line search terminated",
%        "Relative Step Size below tolerance"} |
%        "Deflation operator is Nan/Inf" | "Max Iterations reached" | 
%        "Max Iterations reached" |
%        "Merit line search terminated with rank deficient Jacobian" | 
%        "Deflated Merit line search terminated with rank deficient
%        Jacobian" | "Deflated Merit line search terminated" ]
% Misc Options:
% SysVaryParameters - A structure containing all of the outputs of the
% function [A,A0,scale_x,Ops,SysFixed] = Sys_Input(Sys0,Vary) Used to  speed
%  up the method when calling INS_IEP many times  - [ {[]} | struct ]

if nargin ==4
    Opt = varargin{1};
elseif nargin == 3
    Opt = struct('NDeflations',1);  %Use a default value to avoid error
else
    error('Sys0, Vary and Exp are required inputs')
end
if ~isfield(Exp,"ev")
    error("Please provide the experimental eigenvalues, using Exp.ev")
end


% Process the input of Sys structures, or use precalculated values
if isfield(Opt,'SysVaryParameters')&&isfield(Opt.SysVaryParameters,'A') &&isfield(Opt.SysVaryParameters,'A0')...
        &&isfield(Opt.SysVaryParameters,'x0')&&isfield(Opt.SysVaryParameters,'Ops')&&isfield(Opt.SysVaryParameters,'SysFixed')
    A =        Opt.SysVaryParameters.A;
    A0 =       Opt.SysVaryParameters.A0;
    scale_x =  Opt.SysVaryParameters.x0;
    Ops =      Opt.SysVaryParameters.Ops;
    SysFixed = Opt.SysVaryParameters.SysFixed;
    Opt = rmfield(Opt,'SysVaryParameters');
else
    [A,A0,scale_x,Ops,SysFixed] = Sys_Input(Sys0,varySys);
    if isfield(Opt,'SysVaryParameters')
        Opt = rmfield(Opt,'SysVaryParameters');
    end
end
l = length(scale_x);
%% Set defaults
% INS option defaults:
defaultIEPType = 'Difference';
defaultSysFound = [];
defaultGroundStateFound = [];


%Numerical  option defaults:
% Main method options
defaultMethod ='Good_GN';
defaultNDeflations = 1;
defaultMaxNonMinima = inf;


% Solver options
defaultRegularisation = [];
defaultLinearSolver = "mldivide";
defaultScaled = false;
if length(A{1})<1000
    defaultEigensolver = 'eig';
else
    defaultEigensolver = 'eigs';
end


%Input/output options
defaultVerbose = false;
defaultRecordIterates = true;
defaultSupressConvergedFlag = false;
defaultEigsNotConvergedWarning = true;

%Line search options
defaultLinesearch = 'Armijo';
defaultC1 = 1e-7;
defaultTau = 0.5;
defaultAlpha0 = 1;
defaultMinalpha = 1e-15;
defaultDeflatedLinesearch = false;

%Deflation options
defaultTheta = 2;
defaultSigma = 1;
defaultShiftType = "Multiple";
defaultEpsilon = 0.01;

%Stopping criteria
defaultStepTolerance = 1e-6;
defaultGradientTolerance = 0;
defaultFunctionTolerance = 0;
defaultRelativeStepTolerance = 0;%1e-10;
defaultMaxIter = 1000;
defaultConvergenceFlags = ["Objective less than tolerance","Gradient less than tolerance","Step Size below tolerance","Merit line search terminated","Relative Step Size below tolerance"];

if isfield(Opt,'IEPType')&&Opt.IEPType =="Classic"
    defaultlbylIdentity = speye(length(scale_x)+1);
    isnumerical_lbyl = @(x) isnumeric(x) && all(size(x) == [length(scale_x)+1,length(scale_x)+1]);
else
    defaultlbylIdentity = speye(length(scale_x));
    isnumerical_lbyl = @(x) isnumeric(x) && all(size(x) == [length(scale_x),length(scale_x)]);
end
isSysFoundValid = @(x) isempty(x)||SysCompare(Sys0,varySys,x);

%% Parse Input Parameters

IP = inputParser;
%Validation functions
isnumericscalar = @(x) isscalar(x) && isnumeric(x);
isstringorchar = @(x) isstring(x) || ischar(x);


%Parse required INS variables
addRequired(IP,'Sys0_In')
addRequired(IP,'Vary')
addRequired(IP,'Exp')

%Parse optional INS variables
addParameter(IP,'IEPType',defaultIEPType)
addParameter(IP,'SysFound',defaultSysFound,isSysFoundValid)
addParameter(IP,'GroundStateFound',defaultGroundStateFound)

% Parse Numerical options:
% Main method options
addParameter(IP,'Method',defaultMethod,isstringorchar)
addParameter(IP,'NDeflations',defaultNDeflations,isnumericscalar)
addParameter(IP,'MaxNonMinima',defaultMaxNonMinima,isnumericscalar)

% Solver options
addParameter(IP,'Regularisation',defaultRegularisation)
addParameter(IP,'LinearSolver',defaultLinearSolver,@(x)contains(x,["mldivide","lsqminnorm"]))
addParameter(IP,'Scaled',defaultScaled,@islogical)
addParameter(IP,'Eigensolver',defaultEigensolver)

%Input/output options
addParameter(IP,'Verbose',defaultVerbose,@islogical)
addParameter(IP,'RecordIterates',defaultRecordIterates,@islogical)
addParameter(IP,'SupressConvergedFlag',defaultSupressConvergedFlag,@islogical)
addParameter(IP,'EigsNotConvergedWarning',defaultEigsNotConvergedWarning)

%Line search options
addParameter(IP,'Linesearch',defaultLinesearch,@(str)contains(str,["No","Basic", "Armijo", "Quadratic"]))
addParameter(IP,'c1',defaultC1,isnumericscalar)
addParameter(IP,'tau',defaultTau,isnumericscalar)
addParameter(IP,'alpha0',defaultAlpha0,isnumericscalar)
addParameter(IP,'minalpha',defaultMinalpha,isnumericscalar)
addParameter(IP,'DeflatedLinesearch',defaultDeflatedLinesearch,@islogical)


%Deflation options
addParameter(IP,'theta',defaultTheta,@(x)isnumericscalar(x)||isstringorchar(x))
addParameter(IP,'sigma',defaultSigma,isnumericscalar)
addParameter(IP,'ShiftType',defaultShiftType,@(x)contains(x,["Multiple","Single","Scaled"]))
addParameter(IP,'epsilon',defaultEpsilon,isnumericscalar)
addParameter(IP,'NormWeighting',defaultlbylIdentity,isnumerical_lbyl)

addParameter(IP,'StepTolerance',defaultStepTolerance,isnumericscalar)
addParameter(IP,'GradientTolerance',defaultGradientTolerance,isnumericscalar)
addParameter(IP,'FunctionTolerance',defaultFunctionTolerance,isnumericscalar)
addParameter(IP,'RelativeStepTolerance',defaultRelativeStepTolerance,isnumericscalar)
addParameter(IP,'MaxIter',defaultMaxIter,isnumericscalar)
addParameter(IP,'ConvergenceFlags',defaultConvergenceFlags,isstringorchar)

addParameter(IP,'ScalingMatrix',defaultlbylIdentity,isnumerical_lbyl)
addParameter(IP,'MaxStepSize',[],isnumericscalar)

IP.parse(Sys0,varySys,Exp,varargin{:})
res = IP.Results;

%% Apply INS options

%Turns of the warining from eigs() when it fails to converge for some eigenvalues
if res.EigsNotConvergedWarning
    warning('off','MATLAB:eigs:NotAllEigsConverged')
end

%Checks if Exp.ev is row or column vector
if size(Exp.ev,1)<size(Exp.ev,2)
    constants.ev = Exp.ev';
else
    constants.ev = Exp.ev;
end





% Sets up objective function
if res.IEPType == "Classic"
    obj_fun = @IEP_Evaluate_full;
    if isempty(res.GroundStateFound)
        scale_x(end+1)= -eigs(FormA(scale_x,A,A0),1,'smallestreal');
    else
        scale_x(end+1)= res.GroundStateFound;
    end
    A{end+1} = speye(size(A{1}));
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



% Applies Scaling if used
if res.Scaled
    for i =1:length(scale_x)
        A{i}=scale_x(i)*A{i};
    end
    x0=ones(length(scale_x),1);
else
    x0=scale_x;
end


% Sets constants
constants.A = A;
constants.A0 = A0;
constants.ED = res.Eigensolver;
params.method.constants = constants;


%Parse the input of a previously found systems
if isstruct(res.SysFound)
    %Tests if SysFound was calculated by this package or manually
    if isfield(res.SysFound,'Output')
        [Output(1:length(res.SysFound))] = [res.SysFound(1:length(res.SysFound)).Output];
    else
        %If manually then calculate x and  F(x) for each input system.
        Output = CalculateOutputStructure(res.SysFound,varySys,res.IEPType,constants);
        res.SysFound(1:length(res.SysFound)).Output = Output;
    end
    % If Classic IEP type then check if Groundstate values input
    if res.IEPType == "Classic"%&&~isfield(Output,'GroundStateFound')
        %If not then calculate groundstate values
        for i = 1:length(Output)
            if length(Output(i).DeflatedPoint) == l+1   %groundstate
                if ~isfield(Output,'GroundStateFound')||isempty(Output(i).GroundStateFound)
                    Output(i).GroundStateFound = Output(i).DeflatedPoint(end);
                elseif Output(i).GroundStateFound ~= Output(i).DeflatedPoint(end)
                    error("Mismatch in GroundStateFound and DeflatedPoint in Output")
                end

            elseif isfield(Output,'GroundStateFound')&&~isempty(Output(i).GroundStateFound)
                Output(i).DeflatedPoint(end+1)= Output(i).GroundStateFound;
            elseif isfield(Opt,'GroundStateFound')
                Output(i).DeflatedPoint(end+1) = res.GroundStateFound;
                Output(i).GroundStateFound = res.GroundStateFound;
            else
                smallesteig = -eigs(FormA(Output(i).DeflatedPoint,A(1:end-1),A0),1,'smallestreal');
                Output(i).DeflatedPoint(end+1) = smallesteig;
                Output(i).GroundStateFound = smallesteig;
            end
        end
    elseif res.IEPType == "Difference"
        warn = false;
        for i = 1:length(Output)
            if length(Output(i).DeflatedPoint) == l+1   %groundstate
                if (isfield(Output,'GroundStateFound')&&Output(i).GroundStateFound == Output(i).DeflatedPoint(end))...
                        ||(~isempty(res.GroundStateFound)&& res.GroundStateFound == Output(i).DeflatedPoint(end))
                else
                    warn = true;
                end
                Output(i).DeflatedPoint(end) = [];
            end
        end
        if warn
            warning("SysFound contains solutions with an extra parameter, assuming that previous results calculated using '''Classic''' IEP Type.")
        end
        if ~isempty(res.GroundStateFound)
            params.method.constants.A0 = A0+speye(size(A0))*res.GroundStateFound;
        end
    end
else
    Output = struct("DeflatedPoint",[],"ErrorAtDeflatedPoint",[],"NIter",[],"ConvergenceFlag",[],"FuncCount",[]);
    if res.RecordIterates == true
        Output.Iterates=[];
    end
    if res.IEPType == "Classic"
        Output.GroundStateFound = [];
    end
end




%% Apply numerical method options


%
% Setup Method Parameters
%
params.method.StepMethod = res.Method;
params.method.Verbose = res.Verbose;
params.method.Scaled = res.Scaled;
params.method.Regularisation = res.Regularisation;
params.method.LinearSolver = res.LinearSolver;
params.method.MaxStepSize = res.MaxStepSize;
params.method.ScalingMatrix = res.ScalingMatrix;

if length(res.SysFound) >= res.NDeflations
    error('Please request at least one more minima than input')
end

%Line search Parameters for objective function
params.linesearch.merit.method = res.Linesearch;
params.linesearch.merit.c1 = res.c1;
params.linesearch.merit.tau = res.tau;
params.linesearch.merit.alpha0 = res.alpha0;
params.linesearch.merit.minalpha = res.minalpha;


if params.linesearch.merit.method =="Quadratic"
    if isfield(varargin{1},"tau")
        warning("The optional input Tau has no effect when using the Quadratic linesearch")
    end
end



%Deflation Parameters
params.deflation.theta = res.theta;
params.deflation.sigma = res.sigma;
params.deflation.shifttype = res.ShiftType;
params.deflation.epsilon = res.epsilon;
params.deflation.NormWeighting = res.NormWeighting;

%Convergence Parameters
params.convergence.StepTolerance = res.StepTolerance;
params.convergence.GradientTolerance = res.GradientTolerance;
params.convergence.FunctionTolerance = res.FunctionTolerance;
params.convergence.RelativeStepTolerance = res.RelativeStepTolerance;
params.convergence.MaximumIterations = res.MaxIter;
params.convergence.ConvergenceFlags = res.ConvergenceFlags;

% Set line search objective/merit function and derivatives - phi/gradphi
%First set deflation linesearch
if res.DeflatedLinesearch
    if res.Method ~= "Bad_GN"
        warning("Deflated Linesearch only applied to Bad Gauss-Newton method")
    else
        params.linesearch.Mu = params.linesearch.merit;
        params.linesearch.Mu.phi = @(~,x,constants,DeflatedPts,DeflationParameters) deflation(DeflatedPts,x,DeflationParameters);
        params.linesearch.Mu.gradphi = @(X) X.gradMu;
    end
end
%
% Set objective merit functions and optimisation methods
%
switch res.Method
    case "Newton"
        params.linesearch.merit.phi = @(objective_function,x,constants,DeflatedPts,DeflationParameters) Gradient(objective_function,x,constants,DeflatedPts,DeflationParameters);
        % params.linesearch.merit.phi = @(X) (X.J'*X.R)'*(X.J'*X.R);
        params.linesearch.merit.gradphi = @(X) (X.J'*X.J+X.S)*(X.J'*X.R);
        Fun  = str2func('Deflated_Newton');
    case "Good_GN"
        params.linesearch.merit.phi = @(objective_function,x,constants,~,~) objective_function(x,constants);
        params.linesearch.merit.gradphi = @(X) 2*(X.J'*X.R);
        Fun = str2func('Good_Deflated_GaussNewton');
    case "Bad_GN"
        params.linesearch.merit.phi = @(objective_function,x,constants,~,~) objective_function(x,constants);
        params.linesearch.merit.gradphi = @(X) 2*(X.J'*X.R);
        Fun = str2func('Bad_Deflated_GaussNewton');
    case "RGD_LP"
        Fun  = str2func('RGD_LP');
            if isfield(Opt,"Linesearch")&&Opt.Linesearch ~="No" 
                if isfield(Opt,"Alpha0")&&params.linesearch.merit.alpha0 <=2
                    params.linesearch.merit.method = "No";
                    warning("Line search only applied to Lift and Projection method when Alpha0 > 2.")
                else
                    params.linesearch.merit.phi = @(objective_function,x,constants,~,~) objective_function(x,constants);
                    params.linesearch.merit.gradphi = @(X) 2*(X.J'*X.R);
                end
            else
                params.linesearch.merit.method = "No";
            end
    otherwise
        error('Please select a valid method - Newton/Good_GN/Bad_GN')
end



%% Main Loop

j = 0;

if ~isempty(Output(1).DeflatedPoint)
    if res.Scaled
    DeflatedPts = [Output.DeflatedPoint]./scale_x;
    else
        DeflatedPts = [Output.DeflatedPoint];
    end
else
    DeflatedPts = [];
end
% for  i=length(Output)+notfirstiteration:res.NDeflations
warning('off',"MATLAB:rankDeficientMatrix")
for  i=length(res.SysFound)+1:res.NDeflations
    % for  i=length(Output)+1:res.NDeflations
    if isfield(Output,'DeflatedPoint')&&~isempty([Output.DeflatedPoint])&&any(all(abs([Output.DeflatedPoint] - x0)<1e-16))
        warning('The intial vector is a deflated point, stopping method...')
        break
    end

    % Main optimisation algorithm
    Iterations = []; lastwarn('')
    try
        [Iterations] = Fun(obj_fun,x0,[DeflatedPts,Output(length(res.SysFound)+1:i-1).DeflatedPoint],params,res.RecordIterates);
    catch ME
        msg = getReport(ME);
        warning(msg)
        % stopearly = true;
        break
    end

            [~,ID]=lastwarn;
        if ID == "MATLAB:rankDeficientMatrix"
            warning("Rank deficient matrix detected, consider using the Regularisation option. ")
        end


    if isfield(Output,'GroundStateFound')
        if res.IEPType == "Classic" %Need to fix this option
            Iterations.GroundStateFound = Iterations.DeflatedPoint(end);
            Iterations.GroundStateFound = Iterations.DeflatedPoint(end);
        else
            Iterations.GroundStateFound = [];
        end
    end
    Output(i) = Iterations;

    % Print desired info
    if any(contains(res.ConvergenceFlags, Output(i).ConvergenceFlag))
        j = j+1;
        if ~res.SupressConvergedFlag; OutputNumMinimaFound(j,i-length(res.SysFound));end
    else
        if ~res.SupressConvergedFlag; OutputNumMinimaFound(j,i-length(res.SysFound));end
    end

    % Check for NaNs
    if Output(i).ConvergenceFlag == "NaN/Inf"
        warning("Method diverging to NaN or Inf values, stopped deflations early.")
        break
    elseif Output(i).ConvergenceFlag == "Deflation operator is Nan/Inf"
        warning("Deflation operatoris NaN/Inf, stopped deflations early.")
        break
    elseif (i-length(res.SysFound)-1-j)>res.MaxNonMinima
        warning("Stopping deflations as the maximum number of non minima have been deflated")
        break
    end
end
 
if j ~= i-length(res.SysFound)
    disp("Note: Not all systems found were minima")
end

warning('on',"MATLAB:rankDeficientMatrix")


SysOut1 = [];
if ~isempty(Output(1).DeflatedPoint) && length(Output) - length(res.SysFound) >0
    if res.Scaled
        SysOut1 = Sys_Output([Output(length(res.SysFound)+1:length(Output)).DeflatedPoint].*scale_x,Ops,SysFixed);
    else
        SysOut1 = Sys_Output([Output(length(res.SysFound)+1:length(Output)).DeflatedPoint],Ops,SysFixed);
    end
else
    warning("No new minimising system calculated")
end
warning('on','MATLAB:eigs:NotAllEigsConverged')

for i=1:length(SysOut1)
    SysOut1(i).Output = Output(i+length(res.SysFound));
end
SysOut = [res.SysFound,SysOut1];


end


function f_out = Gradient(objective_function,x,constants,~,~)
[~,R,J] = objective_function(x,constants);
f_out = (J'*R)'*(J'*R);
end

function OutputNumMinimaFound(NumMinFound,OuterIterations)
Deflations = OuterIterations-1;
if Deflations ==1
    if NumMinFound ==1
        disp('Found 1 local minimum in 1 deflation.')
    else
        disp(['Found ',num2str(NumMinFound),' local minima after 1 deflation.'])
    end
else
    if NumMinFound ==1
        disp(['Found 1 local minimum after ',num2str(Deflations),' deflations.'])
    else
        disp(['Found ',num2str(NumMinFound),' local minima after ',num2str(Deflations),' deflations.'])
    end
end
end

function test = SysCompare(Sys0,Vary,SysFound)
test = true;
Varynames = fieldnames(Vary);
for j = length(SysFound)
    for i = 1:length(Varynames)
        if ~all(logical(logical(Vary.(Varynames{i}))+logical(Sys0.(Varynames{i})))==logical(SysFound(j).(Varynames{i})))
            test = fail;
            warning("SysFound structure not compatible")
            return
        end
    end
end
end

function Output = CalculateOutputStructure(SysFound,Vary,IEPType,constants)
for i = 1:length(SysFound)
    [~,~,x]= Sys_Input(SysFound(i),Vary);
    Output(i).DeflatedPoint = x;
    if IEPType == "Classic"
        Output(i).DeflatedPoint(end+1)= -eigs(FormA(Output(i).DeflatedPoint,constants.A(1:end-1),constants.A0),1,'smallestreal');
        Outout(i).GroundStateFound = Output(i).DeflatedPoint(end);
        Output(i).ErrorAtDeflatedPoint = IEP_Evaluate_full(Output(i).DeflatedPoint,constants);
    else
        Output(i).ErrorAtDeflatedPoint = IEP_Evaluate_diff(Output(i).DeflatedPoint,constants);
    end
    Output(i).NIter = 0;
    Output(i).ConvergenceFlag = "";
    Output(i).Iterates = x;
    Output(i).FuncCount = 0;
end
end
