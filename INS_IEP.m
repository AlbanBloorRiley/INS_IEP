function [SysOut] = INS_IEP(Sys0,Vary,Exp,varargin)
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
options = varargin{1};

if length(A{1})<1000
    defaultEigensolver = 'eig';
else
    defaultEigensolver = 'eigs';
end

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

IP.parse(Sys0,Vary,Exp,IEPOptions)
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


IterationFields = ["DeflatedPoint","ErrorAtDeflatedPoint","NIter","ConvergenceFlag","Iterates"];
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
    end

if res.Scaled
    SysOut = Sys_Output([Iterations.DeflatedPoint].*scale_x,Ops,SysFixed);
else
    SysOut = Sys_Output([Iterations.DeflatedPoint],Ops,SysFixed);
end
warning('on','MATLAB:eigs:NotAllEigsConverged')
SysOut = mergestructs(SysOut,Iterations);
end
