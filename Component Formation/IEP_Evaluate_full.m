function [f,varargout] = IEP_Evaluate_full(d,constants)

if ~isfield(constants,'eigenvalueDifferenceSD')
    constants.eigenvalueSD = ones(length(constants.ev),1);
end

if constants.ED =="eig"
    [Q,D] = eig(full(FormA(d,constants.A,constants.A0)),'vector');
    D = D(1:length(constants.ev));
    Q = Q(:,1:length(constants.ev));
else

    [Q,D] = eigs(FormA(d,constants.A,constants.A0),length(constants.ev),'smallestreal');
    D=diag(D);
    if any(isnan(D))
        [~,D1] = eigs(FormA(d,constants.A,constants.A0),max(length(constants.ev)+20,length(constants.A)),'smallestreal');
        D1=diag(D1);
        D = D1(1:length(constants.ev));
    end
    if any(isnan(D))
        [~,D1] = eigs(FormA(d,constants.A,constants.A0),max(2*length(constants.ev)+30,length(constants.A)),'smallestreal');
        D1=diag(D1);
        D = D1(1:length(constants.ev));
    end
    if any(isnan(D))
        warning("eigs() not converged, using eig()");
%         D = eig(FormA(d,constants.A,constants.A0),'vector');
        D = D(1:length(constants.ev));
    end
%     if (D-D1(length(D)))>1e-3
%         disp(D-D1(1:length(D)))
%     end
end
% rcm = 29979.2458;    % reciprocal cm to MHz
% meV = rcm*8.065;
% D = D./meV;
f = sqrt(sum(((D-constants.ev)./constants.eigenvalueSD).^2));
if nargout>1
    varargout{1} = (D-constants.ev)./constants.eigenvalueSD;
end
if nargout>2
    varargout{2} =  FormJ_Lambda(Q,constants.A);
end
if nargout>3
    varargout{3} = FormH_Lambda(Q,constants.A,D,constants.ev);
end
