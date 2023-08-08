function [SysOut, NIter, Flags, Iters, FinalError] = INS_IEP(Sys0,Vary,Exp,varargin)
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

% xint = scale_x;
%  [A1,dint,ev1,Ops1] = SysInput(Sys0);
%     defaultDeflate = "F";

defaultScaled = true;
defaultNMinima = 1;
defaultLinesearch = 'Armijo';
defaultGradientTolerance = 1e-2;
defaultStepTolerance = 1e-8;
defaultObjectiveTolerance = 1e-3;
defaultMaxIter = 10000;
defaultP = 2;
defaultC1 = 1e-4;
defaultAlphaInt = 1;
defaultTau = 0.5;
defaultMinalpha = 1e-4;
defaultTheta = 1;
defaultSingleShift = false;
defaultMethod ='Newton';
defaultDeflateLinesearch = true;
defaultVerbose = false;
defaultIEPType = 'Difference';
defaultSysFound = [];
if length(A{1})<1000
    defaultEigensolver = 'eig';
else
    defaultEigensolver = 'eigs';
end
% defaultA0 = sparse(length(A{1}),length(A{1}));

IP = inputParser;
addRequired(IP,'Sys0_In')
addRequired(IP,'Vary')
addRequired(IP,'Exp')
addParameter(IP,'Method',defaultMethod)

addParameter(IP,'NMinima',defaultNMinima)
addParameter(IP,'Linesearch',defaultLinesearch)

addParameter(IP,'ObjectiveTolerance',defaultObjectiveTolerance)
addParameter(IP,'GradientTolerance',defaultGradientTolerance)
addParameter(IP,'StepTolerance',defaultStepTolerance)

addParameter(IP,'MaxIter',defaultMaxIter)
addParameter(IP,'p',defaultP)
addParameter(IP,'c1',defaultC1)
addParameter(IP,'alphaint',defaultAlphaInt)
addParameter(IP,'tau',defaultTau)
addParameter(IP,'minalpha',defaultMinalpha)
addParameter(IP,'theta',defaultTheta)
addParameter(IP,'SingleShift',defaultSingleShift)
addParameter(IP,'deflatelinesearch',defaultDeflateLinesearch)
addParameter(IP,'Eigensolver',defaultEigensolver)
addParameter(IP,'Verbose',defaultVerbose)
addParameter(IP,'Scaled',defaultScaled)
addParameter(IP,'IEPType',defaultIEPType)

% addParameter(IP,'A0',defaultA0)
addParameter(IP,'SysFound',defaultSysFound)
IP.parse(Sys0,Vary,Exp,varargin{:})
res = IP.Results;

if res.Scaled
    for i =1:length(scale_x)
        A{i}=scale_x(i)*A{i};
    end
    xint=ones(length(scale_x),1);
else
    xint=scale_x;
end



% if ~res.UseInitialGuess
%     xint = ones(length(A),1);
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
    objective_fun = @IEP_Evaluate_full;
    A{end+1} = speye(size(A{1}));
    xint(end+1)=1;
elseif res.IEPType == "Difference"
    objective_fun = @IEP_Evaluate_diff;
else
    error('Please use either the ''Classic'' or ''Difference'' options for IEPType')
end
constants.A = A;
constants.A0 = A0;
constants.ED = res.Eigensolver;




switch res.Method
    case "Newton"
        if res.deflatelinesearch
            phi_0="Mu^2*(Jx'*Rx)'*(Jx'*Rx)";
            phi = @(x,constants,objective_fun,Y,p,theta,SingleShift)F_JDeflated(x,constants,objective_fun,Y,p,theta,SingleShift);
            gradphi = "(Jx'*Rx*gradMu+Mu*(Jx'*Jx+S))*(2*Mu*Jx'*Rx)";
        else
                 phi_0="obj_at_x";
        phi = @(x,constants,objective_fun,Y,p,theta,SingleShift)objective_fun(x,constants);
        gradphi = "2*(Jx'*Rx)";
%             f_0="(Jx'*Rx)'*(Jx'*Rx)";
%             f = @(x,constants,objective_fun,Y,p,theta,SingleShift)objective_fun(x,constants);
%             J_f = "(Jx'*Jx+S)*(Jx'*Rx)";
        end
        Fun= str2func('Deflated_Newton_Method');
    case "Gauss-NewtonT1"
        if res.deflatelinesearch
            if isfield(varargin{1},"deflatelinesearch")
                warning("There is no deflated line search method for the Type 1 Gauss-Newton method. Using an undeflated line search.")
            end
            %             f_0="Mu^2*(Jx'*Rx)'*(Jx'*Rx)";
            %             f = @(x,constants,objective_fun,Y,p,theta,SingleShift)F_JDeflated(x,constants,objective_fun,Y,p,theta,SingleShift);
            %             J_f = "2*(Jx'*Rx*gradMu+Mu*(Jx'*Jx))*(2*Mu*Jx'*Rx)";
        end
        phi_0="obj_at_x";
        phi = @(x,constants,objective_fun,Y,p,theta,SingleShift)objective_fun(x,constants);
        gradphi = "(Jx'*Rx)/dot(Rx,Rx)";

        Fun= str2func('Deflated_GaussNewton_Type1_Method');
    case "Gauss-NewtonT2"
        if res.deflatelinesearch
            phi_0="Mu^2*obj_at_x";
            phi = @(x,constants,objective_fun,Y,p,theta,SingleShift)F_RDeflated(x,constants,objective_fun,Y,p,theta,SingleShift);
            gradphi = "2*Mu*gradMu' + 2*Mu^2*Jx'*Rx";
        else
            
            phi_0="obj_at_x";
            phi = @(x,constants,objective_fun,Y,p,theta,SingleShift)objective_fun(x,constants);
            gradphi = "2*Jx'*Rx";
        end
        Fun= str2func('Deflated_GaussNewton_Type2_Method');
    otherwise
        error('Please select a valid method - Newton/Gauss-NewtonT1/Gauss-NewtonT2')
end
j=0;tic;
for i=length(res.SysFound)+1:res.NMinima
    [Y(:,i),NIter(i),Flags{i},FinalError{i},Iters{i}] = ...
        Fun(objective_fun,xint,constants,res.ObjectiveTolerance,res.GradientTolerance,res.StepTolerance,...
        res.MaxIter,Y,res.p,res.Linesearch,res.c1,res.alphaint,...
        res.tau,res.minalpha,res.theta,res.SingleShift,phi,phi_0,gradphi,res.Verbose);

    if Flags{i} ~="Max Iterations reached"
        j = j+1;
        disp(['Currently found ',num2str(j),' minima in ',num2str(i),' iterations.'])
    else
        disp(['Currently found ',num2str(j),' minima in ',num2str(i),' iterations.'])
    end


    %             xint = Y(:,i)+10*randn(size(Y,1),1);
end
toc


% SysOut = Sys_Output(Y,Ops,Sys0.S);
% SysOut = Sys_Output(Y,Ops,SysFixed);
if res.IEPType == "Classic"
    Y=Y(1:end-1,:);
end
if res.Scaled
    SysOut = Sys_Output(Y.*scale_x,Ops,SysFixed);
else
    SysOut = Sys_Output(Y,Ops,SysFixed);
end


%     for i = length(res.SysFound)+1:res.NMinima
%         Sys(i).Flags = Flags{i};
%         Sys(i).NIter= NIter(i);
%     end
% if ~isempty(res.SysFound)>0
%     Sys(1:length(res.SysFound)) = res.SysFound;
% end
end

function f_out = F_JDeflated(d,constants,objective_fun,Y,p,theta,SingleShift)
M = deflation(Y,d,p,theta,SingleShift);
[~,R,J] = objective_fun(d,constants); f_out = M^2*(J'*R)'*(J'*R);
end
function f_out = F_RDeflated(d,constants,objective_fun,Y,p,theta,SingleShift)
M = deflation(Y,d,p,theta,SingleShift);
f_out = M^2*objective_fun(d,constants);
end
