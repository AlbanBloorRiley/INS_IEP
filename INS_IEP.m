function [SysOut] = INS_IEP(Sys0,Vary,Exp,varargin)
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
    options = varargin{1};
elseif nargin == 3
    options = struct('NDeflations',1);  %Use a default value to avoid error
else  
    error('Sys0, Vary and Exp are required inputs')
end

if ~isfield(Exp,"ev")
    error("Please provide the experimental eigenvalues, using Exp.ev")
end

[A,A0,scale_x,Ops,SysFixed] = Sys_Input(Sys0,Vary);


if length(A{1})<1000
    defaultEigensolver = 'eig';
else
    defaultEigensolver = 'eigs';
end

IEPOptions = [];
IEPOptionFields = ["SysFound","Eigensolver",",EigsNotConvergedWarning","IEPType","Scaled"];
for i = IEPOptionFields
    if isfield(options,i)
        IEPOptions.(i) = options.(i);
        options = rmfield(options,i);
    end
end

mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
% isnumericscalar = @(x)isscalar(x)&&isnumeric(x);
% islbylnumeric = @(x) all(size(x) == [length(x0),length(x0)] )&&isnumeric(x);
% isstringorchar = @(x) isstring(x)||ischar(x);


defaultScaled = false;
defaultIEPType = 'Difference';
defaultSysFound = [];


IP = inputParser;
addRequired(IP,'Sys0_In')
addRequired(IP,'Vary')
addRequired(IP,'Exp')

addParameter(IP,'IEPType',defaultIEPType)
addParameter(IP,'SysFound',defaultSysFound)
addParameter(IP,'Eigensolver',defaultEigensolver)
addParameter(IP,'Scaled',defaultScaled,@islogical)
addParameter(IP,'EigsNotConvergedWarning',true)

if isempty(IEPOptions)
    IP.parse(Sys0,Vary,Exp)
else
    IP.parse(Sys0,Vary,Exp,IEPOptions)
end
res = IP.Results;

%Turns of the warining from eigs() when it fails to converge for some eigenvalues
if res.EigsNotConvergedWarning
    warning('off','MATLAB:eigs:NotAllEigsConverged')
end



if res.Scaled
    for i =1:length(scale_x)
        A{i}=scale_x(i)*A{i};
    end
    x0=ones(length(scale_x),1);
else
    x0=scale_x;
end


IterationFields = ["DeflatedPoint","ErrorAtDeflatedPoint","NIter","ConvergenceFlag","Iterates","GroundStateFound"];
if isstruct(res.SysFound)
    for j = IterationFields
        if isfield(res.SysFound,j)
            PreviouslyFoundIterations.(j) = res.SysFound.(j);
            res.SysFound = rmfield(res.SysFound,j);
        end
    end
    options.PreviouslyFoundIterations = PreviouslyFoundIterations;

end

%Checks if Exp.ev is row or column vector
if size(Exp.ev,1)<size(Exp.ev,2)
    constants.ev = Exp.ev';
else
    constants.ev = Exp.ev;
end


if res.IEPType == "Classic"
    obj_fun = @IEP_Evaluate_full;
    A{end+1} = speye(size(A{1}));
    x0(end+1)=1;
    if isstruct(res.SysFound)%
       if ~isfield(PreviouslyFoundIterations,"GroundStateFound")
           error("Please provide the value of the ground state found for previously found systems")
       end
        for i = 1:length(PreviouslyFoundIterations)
           PreviouslyFoundIterations(i).DeflatedPoint(end+1)=PreviouslyFoundIterations(i).GroundStateFound;
       end
    end
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
options.constants = constants;


Iterations = DMin(obj_fun,x0,options);


if res.IEPType == "Classic" %Need to fix this option
    for i = 1:length(Iterations)
        GroundStateFound(i).GroundStateFound = Iterations(i).DeflatedPoint(end);
        Iterations(i).DeflatedPoint(end) = [];
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
end
