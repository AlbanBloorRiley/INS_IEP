function [SysOut,options,params,obj_fun,Iterations] = INS_IEP(Sys0,Vary,Exp,varargin)
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
% Vary should contain all the same fields of Sys (not including S), and any
% entry that is non-zero specifies that it is a parameter to be varied,
% using the initial value given in Sys. Note that parameters that are given
% by Sys but are not to be varied will still be used, just as constants.
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
%OPTIONAL PARAMETERS
%IEP:
%IEPType - Defines which formulation of the IEP to use - [ Classic |
%       {Difference} ]
%Eigensolver - Specifies which eigensolver to use - [ eig | eigs ] (default
%       is decided based on the size of the matrices used)
%Scaled - specifies if the parameters are scaled to the size of the initil
%       values - [ true | {false} ]
%EigsNotConvergedWarning - Specifies if a warning should be output of the
%       eigensolver fails to converge for some of the eigenvalues - [ {true} |
%       false ]
%NUMERICAL:
%NDeflations - The number of local minima you wish to find, the output -
% [ {1} | positive integer ]
%Method - The optimisation method desired - [ {Newton} | GaussNewtonT1 |
%           GaussNewtonT2 ]
%MaxIter - The integer value  of Maximum iterations per deflation -
%           [ {1000} | positive integer ]
%Linesearch - line search method - [ No | {Armijo} | Quadratic ]
%theta - The deflation exponent - [ {2} | positive integer | exp ]
%c1 - the armijo line search parameter - [ {1e-4} | positive scalar ]
%alpha0 - Initial value of the line search parameter each iteration - [ 1
%           | positive scalar ]
%Tau - The value of the decrease in the line search parameter - [ {0.5} |
%         scalar in [0,1] ]
%Minalpha - The minimum value of the line search parameter - [ 1e-18 |
%       positive scalar ]
%Verbose - Output value of objective function and gradient at each
%            iteration - [ true | {false} ]
%Constants - Any constants the objective funtion provided requires - [ {}
%        | struct ]
%PreviouslyFoundIterations - A struct of prevously found points to be deflated
%        - [ {} | l \times x matrix ], where x is the number of previous
%        deflations
%ConvergenceFlags - The flags that are considered to mean convergence - [
%        {Objective less than tolerance} | {Gradient less than tolerance} |
%       Step Size too small | Line search failed | Max Iterations reached |
%        Divergence Detected | NaN ]
%Sigma - The value of the shift of the deflation - [ {1} | scalar ]
%SingleShift - Deflation shift should only be applied once [ true |
%        {false} ]
%Epsilon - Tolerence on the application of deflation or line search - [
%        {0.01} | scalar in [0,1]]
%NormWeighting - Weighting applied to the norms in calculation of
%           deflation operators - [ {} | l \times l matrix ]
%GradientTolerance - Stopping criterion based on gradient of function -
%        [ {0} | positive scalar ]
%StepTolerance - Stopping criterion based on difference of consecutive
%        steps - [ {1e-8} | positive scalar ]
%ObjectiveTolerance - Stopping criterion based on value of objective
%        function - [ {0} | positive scalar ]
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




%Precalculated INS operators
if isfield(Opt,'INSOperators')&&isfield(Opt.INSOperators,'A') &&isfield(Opt.INSOperators,'A0')...
        &&isfield(Opt.INSOperators,'x0')&&isfield(Opt.INSOperators,'Ops')&&isfield(Opt.INSOperators,'SysFixed')
    A =        Opt.INSOperators.A;
    A0 =       Opt.INSOperators.A0;
    scale_x =  Opt.INSOperators.x0;
    Ops =      Opt.INSOperators.Ops;
    SysFixed = Opt.INSOperators.SysFixed;
    Opt = rmfield(Opt,'INSOperators');
else
    [A,A0,scale_x,Ops,SysFixed] = Sys_Input(Sys0,Vary);
    if isfield(Opt,'INSOperators')
        Opt = rmfield(Opt,'INSOperators');
    end
end
if length(A{1})<1000
    defaultEigensolver = 'eig';
else
    defaultEigensolver = 'eigs';
end
% IEPOptions = [];
% IEPOptionFields = ["SysFound","Eigensolver",",EigsNotConvergedWarning","IEPType","Scaled","GroundStateFound"];
% for i = IEPOptionFields
%     if isfield(Opt,i)
%         IEPOptions.(i) = Opt.(i);
%         Opt = rmfield(Opt,i);
%     end
% end

% INS option defaults
defaultScaled = false;
defaultIEPType = 'Difference';
defaultSysFound = [];

%Numerical Method  option defaults
defaultMethod ='Good_GN';
defaultNDeflations = 1;
defaultVerbose = false;
defaultSupressConvergedFlag = false;
% defaultConstants=[];
defaultRegularisation = [];
defaultLinearSolver = "mldivide";

defaultPreviouslyFoundIterations = [];
defaultConvergenceFlags = ["Objective less than tolerance","Gradient less than tolerance","Step Size below tolerance","Merit line search terminated"];
defaultRecordIterates = true;
defaultMaxNonMinima = 10;
defaultRecordTimes = false;

defaultLinesearch = 'Armijo';
defaultC1 = 1e-4;
defaultTau = 0.5;
defaultAlpha0 = 1;
defaultMinalpha = 1e-18;
defaultMuLinesearch = "No";
defaultDeflatedLinesearch = "No";

defaultTheta = 2;
defaultSigma = 1;
defaultShiftType = "Multiple";
defaultEpsilon = 0.01;
defaultNormWeighting = speye(length(scale_x));

defaultStepTolerance = 1e-6;
defaultGradientTolerance = 0;
defaultObjectiveTolerance = 0;
defaultMaxIter = 1000;

if isfield(Opt,'IEPType')&&Opt.IEPType =="Classic"
    defaultlbylIdentity = speye(length(scale_x)+1);
    isnumerical_lbyl = @(x) isnumeric(x) && all(size(x) == [length(scale_x)+1,length(scale_x)+1]);
else
    defaultlbylIdentity = speye(length(scale_x));
    isnumerical_lbyl = @(x) isnumeric(x) && all(size(x) == [length(scale_x),length(scale_x)]);
end

%
% Parse Input Parameters
%

%Validation functions
mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
% islbylnumeric = @(x) all(size(x) == [length(x0),length(x0)] )&&isnumeric(x);

isnumericscalar = @(x) isscalar(x) && isnumeric(x);
isstringorchar = @(x) isstring(x) || ischar(x);
isvalidstruct =@(x) (isstruct(x)&&4==sum(contains(fieldnames(x),[{'DeflatedPoint'}, ...
    {'ErrorAtDeflatedPoint'},{'NIter'},{'ConvergenceFlag'}])))||isempty(x);



IP = inputParser;

%Add required INS variables
addRequired(IP,'Sys0_In')
addRequired(IP,'Vary')
addRequired(IP,'Exp')

%Add optional INS variables
addParameter(IP,'IEPType',defaultIEPType)
addParameter(IP,'SysFound',defaultSysFound)
addParameter(IP,'Eigensolver',defaultEigensolver)
addParameter(IP,'Scaled',defaultScaled,@islogical)
addParameter(IP,'EigsNotConvergedWarning',true)
addParameter(IP,'GroundStateFound',[])
addParameter(IP,'OrthogonaliseBasis',true,@islogical)

%Add optional Numerical method variables

addParameter(IP,'Method',defaultMethod,isstringorchar)
addParameter(IP,'NDeflations',defaultNDeflations,isnumericscalar)
addParameter(IP,'Verbose',defaultVerbose,@islogical)
addParameter(IP,'SupressConvergedFlag',defaultSupressConvergedFlag,@islogical)
% addParameter(IP,'constants',defaultConstants)
addParameter(IP,'Regularisation',defaultRegularisation)
addParameter(IP,'LinearSolver',defaultLinearSolver,@(x)contains(x,["mldivide","lsqminnorm"]))

addParameter(IP,'PreviouslyFoundIterations',defaultPreviouslyFoundIterations,isvalidstruct)
addParameter(IP,'ConvergenceFlags',defaultConvergenceFlags,isstringorchar)
addParameter(IP,'RecordIterates',defaultRecordIterates,@islogical)
addParameter(IP,'MaxNonMinima',defaultMaxNonMinima,isnumericscalar)
addParameter(IP,'RecordTimes',defaultRecordTimes,@islogical)

addParameter(IP,'Linesearch',defaultLinesearch,isstringorchar)
addParameter(IP,'c1',defaultC1,isnumericscalar)
addParameter(IP,'tau',defaultTau,isnumericscalar)
addParameter(IP,'alpha0',defaultAlpha0,isnumericscalar)
addParameter(IP,'minalpha',defaultMinalpha,isnumericscalar)
addParameter(IP,'MuLinesearch',defaultMuLinesearch,isstringorchar)
addParameter(IP,'Mualpha0',defaultAlpha0,isnumericscalar)

addParameter(IP,'DeflatedLinesearch',defaultDeflatedLinesearch,isstringorchar)
addParameter(IP,'Deflatedalpha0',[],isnumericscalar)
addParameter(IP,'Deflatedc1',[],isnumericscalar)
addParameter(IP,'Deflatedminalpha',[],isnumericscalar)
addParameter(IP,'Deflatedtau',[],isnumericscalar)

addParameter(IP,'theta',defaultTheta,@(x)isnumericscalar(x)||isstringorchar(x))
addParameter(IP,'sigma',defaultSigma,isnumericscalar)
addParameter(IP,'ShiftType',defaultShiftType,@(x)contains(x,["Multiple","Single","Scaled"]))
addParameter(IP,'epsilon',defaultEpsilon,isnumericscalar)
addParameter(IP,'NormWeighting',defaultlbylIdentity,isnumerical_lbyl)

addParameter(IP,'ObjectiveTolerance',defaultObjectiveTolerance,isnumericscalar)
addParameter(IP,'GradientTolerance',defaultGradientTolerance,isnumericscalar)
addParameter(IP,'StepTolerance',defaultStepTolerance,isnumericscalar)
addParameter(IP,'MaxIter',defaultMaxIter,isnumericscalar)

addParameter(IP,'ScalingMatrix',defaultlbylIdentity,isnumerical_lbyl)
addParameter(IP,'MaxStepSize',[],isnumericscalar)

% if isempty(IEPOptions)
%     IP.parse(Sys0,Vary,Exp)
% else
IP.parse(Sys0,Vary,Exp,varargin{:})
% end
res = IP.Results;


%% Apply INS options

%Turns of the warining from eigs() when it fails to converge for some eigenvalues
if res.EigsNotConvergedWarning
    warning('off','MATLAB:eigs:NotAllEigsConverged')
end

% Parses previously found results
IterationFields = ["DeflatedPoint","ErrorAtDeflatedPoint","NIter","ConvergenceFlag","Iterates"];
if isstruct(res.SysFound)
    for j = IterationFields
        if isfield(res.SysFound,j)
            for i = 1:length(res.SysFound)
                Iterations(i).(j) = res.SysFound(i).(j);
            end
            res.SysFound = rmfield(res.SysFound,j);
        end
    end

    if isfield(res.SysFound,"GroundStateFound")&&res.IEPType == "Classic"
        if isfield(res.SysFound,"GroundStateFound")
            Iterations.DeflatedPoint(end+1) =  res.SysFound.GroundStateFound;
            % res.SysFound = rmfield(res.SysFound,'GroundStateFound')
        elseif ~isempty(res.GroundStateFound)
            Iterations.DeflatedPoint(end+1) =  res.GroundStateFound;
        else
            error("Ground state of previous SysFound must be provided")
        end
    end
    notfirstiteration = 1;
    % Opt.PreviouslyFoundIterations = Iterations;
else
    notfirstiteration = 0;
    Iterations = struct("DeflatedPoint",[],"ErrorAtDeflatedPoint",[],"NIter",[],"ConvergenceFlag",[],"FuncCount",[]);
    if res.RecordIterates == true
        Iterations.Iterates=[];
    end
    % Iterations = [];
end

%Checks if Exp.ev is row or column vector
if size(Exp.ev,1)<size(Exp.ev,2)
    constants.ev = Exp.ev';
else
    constants.ev = Exp.ev;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to do this

% if isfield(Opt,"Method")&& Opt.Method == "LP"
%     if res.IEPType =="Difference"
%         warning("If using the Lift and Project method then the Classic IEP Type will be used")
%         res.IEPType = "Classic";
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sets up objective function
if res.IEPType == "Classic"
    obj_fun = @IEP_Evaluate_full;

    % scale_x(end+1)=1;
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


% Applies Scaling
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
% Opt.constants = constants;
params.method.constants = constants;


%
% if isfield(Opt,"Method")&& Opt.Method == "LP"
%     Opt.Method = "GradientDescent";
%     B = zeros(length(A));
%     for i = 1:(length(A))
%         for j = 1:(length(A))
%             B(i,j) = sum(sum(A{j}'.*A{i}));
%         end
%     end
%     if isfield(Opt,"ScalingMatrix")
%         Opt.ScalingMatrix = Opt.ScalingMatrix*2*inv(B);
%     else
%         Opt.ScalingMatrix = 2*inv(B);
%     end
% end





%% Apply numerical method options

%
% Setup Method Parameters
%
params.method.StepMethod = res.Method;
params.method.Verbose = res.Verbose;
params.method.Scaled = res.Scaled;
% params.method.constants = res.constants;
params.method.Regularisation = res.Regularisation;
params.method.LinearSolver = res.LinearSolver;
params.method.MaxStepSize = res.MaxStepSize;
params.method.ScalingMatrix = res.ScalingMatrix;

if length(res.PreviouslyFoundIterations) >= res.NDeflations
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


%Line search Parameters for deflated objective function
if res.Method=="Bad_GN"||res.Method=="GradientDescent"||res.Method=="LP"
    if isfield(varargin{1},"LinearSolver")
        warning("The optional input LinearSolver has no effect when using the ""Bad"" Deflated Gauss Newton or deflated Gradient Descent method");
    end
    if isempty(res.DeflatedLinesearch)
        params.linesearch.deflatedmerit.method = res.Linesearch;
    else
        params.linesearch.deflatedmerit.method = res.DeflatedLinesearch;
        if params.linesearch.deflatedmerit.method =="Quadratic"
            if ~isempty(res.Deflatedtau)
                warning("The optional input DeflatedTau has no effect when using the deflated Quadratic linesearch")
            end
        end
    end
    for i = ["c1","tau","alpha0","minalpha"]
        if isempty(res.(strcat("Deflated",i)))
            params.linesearch.deflatedmerit.(i) = res.(i);
        else
            params.linesearch.deflatedmerit.(i) = res.(strcat("Deflated",i));
        end
    end
elseif (res.DeflatedLinesearch)~="No"
    warning("Deflated linesearch is only applied to the ""Bad"" deflated Gauss-Newton method")
end

%parameters for line search on Mu, not all adaptable currently
params.linesearch.Mu.method = res.MuLinesearch;
params.linesearch.Mu.c1 = res.c1;
params.linesearch.Mu.tau = res.tau;
params.linesearch.Mu.alpha0 = res.Mualpha0;
params.linesearch.Mu.minalpha = res.minalpha;



%Deflation Parameters
params.deflation.theta = res.theta;
params.deflation.sigma = res.sigma;
params.deflation.shifttype = res.ShiftType;
params.deflation.epsilon = res.epsilon;
params.deflation.NormWeighting = res.NormWeighting;

%Convergence Parameters
params.convergence.ObjectiveTolerance = res.ObjectiveTolerance;
params.convergence.GradientTolerance = res.GradientTolerance;
params.convergence.StepTolerance = res.StepTolerance;
params.convergence.MaximumIterations = res.MaxIter;


% Set line search objective/merit function and derivatives - phi/gradphi
%First set deflation linesearch
params.linesearch.Mu.phi = @(~,x,constants,DeflatedPts,DeflationParameters) deflation(DeflatedPts,x,DeflationParameters);
params.linesearch.Mu.gradphi = @(X) X.gradMu;

%
% Set objective merit functions and optimisation methods
%
switch res.Method
    case "Newton"
        params.linesearch.merit.phi = @(objective_function,x,constants,DeflatedPts,DeflationParameters) Gradient(objective_function,x,constants,DeflatedPts,DeflationParameters);
        params.linesearch.merit.gradphi = @(X) (X.J'*X.J+X.S)*(X.J'*X.R);
        if isfield(params.linesearch,"deflatedmerit")&&params.linesearch.deflatedmerit.method ~= "No"
            params.linesearch.deflatedmerit.phi = @(objective_function,x,constants,DeflatedPts,DeflationParameters) Deflated_Gradient(objective_function,x,constants,DeflatedPts,DeflationParameters);
            params.linesearch.deflatedmerit.gradphi = @(X) 2*X.Mu*X.gradMu'*(X.J'*X.R)'*(X.J'*X.R)+ 2*X.Mu^2*(X.J'*X.J+X.S)*(X.J'*X.R);
        end
        Fun  = str2func('Deflated_Newton');
    case "Good_GN"
        if isfield(params.linesearch,"deflatedmerit")&&params.linesearch.deflatedmerit.method ~= "No"
            warning("Note there is no deflated line search for any method other than for the '''Bad''' deflated Gauss-Newton and '''Newton''' methods.")
        end
        params.linesearch.merit.phi = @(objective_function,x,constants,~,~) objective_function(x,constants);
        params.linesearch.merit.gradphi = @(X) 2*(X.J'*X.R);
        Fun = str2func('Good_Deflated_GaussNewton');
    case "Bad_GN"
        params.linesearch.merit.phi = @(objective_function,x,constants,~,~) objective_function(x,constants);
        params.linesearch.merit.gradphi = @(X) 2*(X.J'*X.R);
        if params.linesearch.deflatedmerit.method ~= "No"
            params.linesearch.deflatedmerit.phi = @(objective_function,x,constants,DeflatedPts,DeflationParameters) deflation(DeflatedPts,x,DeflationParameters)^2*objective_function(x,constants);
            params.linesearch.deflatedmerit.gradphi = @(X) 2*X.Mu^2*(X.J'*X.R)+2*X.Mu*X.R'*X.R;
        end
        Fun = str2func('Bad_Deflated_GaussNewton');
    case "LP"
        Fun  = str2func('LP');
        if res.NDeflations >1 || isstruct(res.SysFound)
            warning("The Lift and Project method does not curently suport deflation, only one minimum can me requested")
            res.NDeflations=1;
            Opt.PreviouslyFoundIterations = [];
        else
            if isfield(Opt,"Linesearch")&&params.linesearch.merit.method ~= "No"||isfield(params.linesearch,"deflatedmerit")&&params.linesearch.deflatedmerit.method ~= "No"
                warning("Note that no line search can be used with the Lift and Projection method.")
            end

            if res.IEPType ~= "Classic"
                warning("When using the LP method the ""Classic"" IEP type must be used, switching to the  RGD_LP method.")
                res.Method = "RGD_LP";
                Fun  = str2func('RGD_LP');
                if res.NDeflations >1 || isstruct(res.SysFound)
                    warning("The Lift and Project method does not curently suport deflation")
                    res.NDeflations=1;
                    Opt.PreviouslyFoundIterations = [];
                    % Opt = rmfield(Opt, "PreviouslyFoundIterations");
                end

            end
        end
    case "RGD_LP"
        Fun  = str2func('RGD_LP');
        if res.NDeflations >1 || isstruct(res.SysFound)
            warning("The Lift and Project method does not curently suport deflation, only one minimum can me requested") 
            res.NDeflations=1;
            Opt.PreviouslyFoundIterations = [];
        else
            if isfield(Opt,"Linesearch")&&params.linesearch.merit.method ~= "No"||isfield(params.linesearch,"deflatedmerit")&&params.linesearch.deflatedmerit.method ~= "No"
                warning("Note that no line search can be used with the Lift and Projection method.")
            end

            if res.NDeflations >1 || isstruct(res.SysFound)
                warning("The Lift and Project method does not curently suport deflation, only one minimum can me requested")
                res.NDeflations=1;
                Opt.PreviouslyFoundIterations = [];
            end
        end
    otherwise
        error('Please select a valid method - Newton/Good_GN/Bad_GN')
end



%% Main Loop

j = 0;
if res.RecordTimes
    % Times = zeros([res.NDeflations,1]);
    % stopearly=false;
    Times{1,:,:} = [];
    tic
end

for  i=length(Iterations)+notfirstiteration:res.NDeflations
    % for  i=length(Iterations)+1:res.NDeflations
    if isfield(Iterations,'DeflatedPoint')&&~isempty([Iterations.DeflatedPoint])&&any(all(abs([Iterations.DeflatedPoint] - x0)<1e-16))
        warning('The intial vector is a deflated point, stopping method...')
        break
    end

    % Main optimisation algorithm
    try
        [Iterations(i)] = Fun(obj_fun,x0,[Iterations.DeflatedPoint],params,res.RecordIterates);
    catch ME
        msg = getReport(ME);
        warning(msg)
        % stopearly = true;
        break
    end
    %Record Time to compute if requested
    if res.RecordTimes
        Times{i,:,:} = toc;
        tic
    end

    % Print desired info
    if ~res.SupressConvergedFlag
        if any(contains(res.ConvergenceFlags, Iterations(i).ConvergenceFlag))
            j = j+1;
            OutputNumMinimaFound(j,i-length(res.PreviouslyFoundIterations))
        else
            OutputNumMinimaFound(j,i-length(res.PreviouslyFoundIterations))
        end
    end
    % Check for NaNs
    if Iterations(i).ConvergenceFlag == "NaN/Inf"
        warning("Method diverging to NaN or Inf values, stopped deflations early.")
        % stopearly = true;
        break
    end
    if Iterations(i).ConvergenceFlag == "Deflation operator is Nan/Inf"
        warning("Deflation operatoris NaN/Inf, stopped deflations early.")
        % stopearly = true;
        break
    end
    if (i-length(res.PreviouslyFoundIterations)-1-j)>res.MaxNonMinima
        warning("Stopping deflations as the maximum number of non minima have been deflated")
        % stopearly = true;
        break
    end
end
if res.RecordTimes
    % if stopearly
    %     Times{i,:,:} = NaN;
    % end
    Iterations = cell2struct([struct2cell(Iterations);reshape(Times,1,1,length(Iterations))],[fieldnames(Iterations);'Times']);
end















if res.IEPType == "Classic" %Need to fix this option
    for i = 1:length(Iterations)
        if res.Scaled
            GroundStateFound(i).GroundStateFound = Iterations(i).DeflatedPoint(end)*scale_x(end);
        else
            GroundStateFound(i).GroundStateFound = Iterations(i).DeflatedPoint(end);
        end
        Iterations(i).DeflatedPoint(end) = [];
    end
    if res.Scaled
        scale_x(end)=[];
    end
end
if res.Scaled
    SysOut = Sys_Output([Iterations.DeflatedPoint].*scale_x,Ops,SysFixed);
else
    SysOut = Sys_Output([Iterations.DeflatedPoint],Ops,SysFixed);
end
warning('on','MATLAB:eigs:NotAllEigsConverged')
SysOut = mergestructs(SysOut,Iterations);
if res.IEPType == "Classic"
    SysOut = mergestructs(SysOut,GroundStateFound);
end

options.Eigensolver = res.Eigensolver;
options.EigsNotConvergedWarning = res.EigsNotConvergedWarning;
options.Exp = res.Exp;
options.IEPType = res.IEPType;
options.Scaled = res.Scaled;
options.Sys0_In = res.Sys0_In;
options.SysFound = res.SysFound;
options.Vary = res.Vary;
end


% function [Iterations,res,params] = DMin(obj_fun,x0,varargin)
% %DMIN Finds multiple minima of the problem using deflation.
% %
% % Y= DMIN(obj_fun,x0)   Outputs a structure that contains the points that
% % minimise the problem 'obj_fun' initialised as x0and the number of Iterations till
% % convergence; the type of 'convergence' achieved; all the iterations
% % made until convergence; the finall error at convergence. Where
% % obj_fun(x,constants) is a function handle that returns the current
% %    error of the least squares formulation, given the current guess, x,
% %    and optional parameters, constants. It must also  output the
% %    Residual and the Jacbian and (optionally) Hessian (if using the Newton
% % method) of said residual of the least squares problem: [F,R,J,H] =
% % obj_fun(x,constants)
% %
% %
% %DMIN PARAMETERS
% %NDeflations - The number of local minima you wish to find, the output -
% % [ {1} | positive integer ]
% %Method - The optimisation method desired - [ {Newton} | GaussNewtonT1 |
% %           GaussNewtonT2 ]
% %MaxIter - The integer value  of Maximum iterations per deflation -
% %           [ {1000} | positive integer ]
% %Linesearch - line search method - [ No | Basic | {Armijo} | Quadratic ]
% %theta - The deflation exponent - [ {2} | positive integer | exp ]
% %c1 - the armijo line search parameter - [ {1e-4} | positive scalar ]
% %alpha0 - Initial value of the line search parameter each iteration - [ 1
% %           | positive scalar ]
% %Tau - The value of the decrease in the line search parameter - [ {0.5} |
% %         scalar in [0,1] ]
% %Minalpha - The minimum value of the line search parameter - [ 1e-10 |
% %       positive scalar ]
% %Verbose - Output value of objective function and gradient at each
% %            iteration - [ true | {false} ]
% %Scaled - Scale the variables by the value of the initial guess - [ true |
% %       {false} ]
% %Constants - Any constants the objective funtion provided requires - [ {}
% %        | struct ]
% %Regularisation - The value regularing parameter - [ {0} ]
% %LinearSolver - The linear solver used - [ {mldivide} | lsqminnorm ]
% %PreviouslyFoundIterations - A struct of prevously found points to be deflated
% %        - [ {} | l \times x matrix ], where x is the number of previous
% %        deflations
% %ConvergenceFlags - The flags that are considered to mean convergence - [
% %        {Objective less than tolerance} | {Gradient less than tolerance} |
% %       Step Size too small | Line search failed | Max Iterations reached |
% %        Divergence Detected | NaN ]
% %DeflateLinesearch - Line searches should be applied to the deflated
% %        system - [ {true} | false ]
% %Sigma - The value of the shift of the deflation - [ {1} | scalar ]
% %SingleShift - Deflation shift should only be applied once [ true |
% %        {false} ]
% %Epsilon - Tolerence on the application of deflation or line search - [
% %        {0.01} | scalar ]
% %NormWeighting - Weighting applied to the norms in calculation of
% %           deflation operators - [ {} | l \times l matrix ]
% %GradientTolerance - Stopping criterion based on gradient of function -
% %        [ {0} | positive scalar ]
% %StepTolerance - Stopping criterion based on difference of consecutive
% %        steps - [ {1e-8} | positive scalar ]
% %ObjectiveTolerance - Stopping criterion based on value of objective
% %        function - [ {0} | positive scalar ]
%
% Auxiliary functions

function f_out = Deflated_Gradient(objective_function,x,constants,DeflatedPts,DeflationParameters)
[~,R,J] = objective_function(x,constants);
f_out = deflation(DeflatedPts,x,DeflationParameters)^2*(J'*R)'*(J'*R);
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
