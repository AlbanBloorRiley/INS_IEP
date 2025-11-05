function [f,varargout] = IEP_Evaluate_Zerodiff(x,constants)
% IEP_Evaluate Evaluates the inverse eigenvalue (difference) problem.
%f = IEP_Evaluate(x,constants) calculates the error between the current and
%prescribed eigenvalue (differences).
%    * x is the current parameter guess where A = A_0 + sum(A_i*x_i)
%    * constants is a structure containing the values of the constants A,
%      A0, the type of eigendecompostion to be used and the prescribed 
%      eigenvalues. Stored as constants.A,
%      constants.A0, constants.ED ,constants.ev (a column vector), and 
%      optionally the standard deviation of the eigenvalue differences.
%
%   [f,J] = IEP_Evaluate(x,constants) also outputs the residual function                      
%
%   [f,R,J] = IEP_Evaluate(x,constants) also outputs the Jacobian of R
%   
%   [f,R,J,H] = IEP_Evaluate(x,constants) also outputs the Hessian of R.                     
if ~isfield(constants,'eigenvalueDifferenceSD')
    constants.eigenvalueDifferenceSD = ones(length(constants.ev)-1,1);
end
if isfield(constants,'ED') && constants.ED =="eig"
    [Q,D] = eig(full(FormA(x,constants.A,constants.A0)),'vector');
    D = D(1:length(constants.ev));
    Q = Q(:,1:length(constants.ev));
else

    [Q,D1] = eigs(FormA(x,constants.A,constants.A0),length(constants.ev),'smallestreal');
    D1=diag(D1);
        D = (D1(1:length(constants.ev)));
    % if any(isnan(D))
    %     [Q,D1] = eigs(FormA(x,constants.A,constants.A0),2*length(constants.ev),'smallestreal','subspacedimension',2*min([length(constants.ev)+10,2*length(constants.ev),length(constants.A{1})]));
    %     D1=diag(D1);
    %     D = (D1(1:length(constants.ev)));
    % end
    if any(isnan(D))
        [Q,D1] = eigs(FormA(x,constants.A,constants.A0),3*length(constants.ev),'smallestreal','subspacedimension',2*min([max(length(constants.ev)+50,3*length(constants.ev)),length(constants.A{1})]));
        D1=diag(D1);
        D = D1(1:length(constants.ev));
    end
    if any(isnan(D))
        warning("eigs() not converged, using eig()");
        [Q,D] = eig(full(FormA(x,constants.A,constants.A0)),'vector');
        D = D(1:length(constants.ev));
    end
%     if (D-D1(length(D)))>1e-3
%         disp(D-D1(1:length(D)))
%     end
end
if isfield(constants,'MintStruct')
    constants.MintStruct.Opt.Eigs = D;
    constants.MintStruct.Opt.Vecs = Q;
    cross_sect = mint(constants.MintStruct.Sys,constants.MintStruct.Exp,constants.MintStruct.Opt);
    plot(constants.MintStruct.Exp.Energy,cross_sect*1e-4,'linewidth',1.2)
end

% R = (D(2:end)-D(1:end-1))-(constants.ev(2:end)-constants.ev(1:end-1));
% J = FormJ_Lambda(Q(:,2:end),constants.A) - FormJ_Lambda(Q(:,1:end-1),constants.A);
%  L = D-constants.ev;
    RL = ((D(2:end) - D(1)) - (constants.ev(2:end) - constants.ev(1)))./constants.eigenvalueDifferenceSD;
    f = sqrt(sum((RL).^2));
if nargout>1
    varargout{1} = RL;
end
if nargout>2
    LJ = FormJ_Lambda(Q(:,1:length(constants.ev)),constants.A);
%     LJ = FormJ_Lambda_old_1(Q(:,1:length(constants.ev)),constants.A);

    varargout{2} = LJ(2:end,:) - LJ(1,:);
end
if nargout>3
    LH =  FormH_Lambda(Q,constants.A,D,constants.ev);
    varargout{3} = LH(:,:,2:end) - LH(:,:,1);
end
