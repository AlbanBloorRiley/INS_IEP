function [Mu,gradMu] = deflation(Y,x,varargin)
%    DEFLATION Calculates shifted deflation operators
%    Mu = DEFLATION(y,x) produces the shifted deflation operator Mu at point
%    x, given a row vector of deflated points Y. 
%
%    mu(x) = prod_i (i
%    The default operator uses the
%    square of the 2-norm and a shift of 1.
%
%    [Mu,gradMu] = DEFLATION(y,x) also produces the gradient of the deflation
%    operator.
%
%    [Mu,gradMu] = DEFLATION(y,x,theta) changes the exponent on the norm to
%    theta. Unless theta is 'exp', then it calculates the exponential deflation
%    operator.
%
%    [Mu,gradMu] = DEFLATION(y,x,theta,sigma) Changes the value of the shift to
%    sigma
%
%    [Mu,gradMu] = DEFLATION(y,x,theta,sigma,true) Uses the single shift
%    strategy for multiple deflations.
%
%    [Mu,gradMu] = DEFLATION(y,x,DeflationParameters) An alternative way to
%    pass in the optional variables, given that all options are included.

if isempty(Y)
    Mu = 1;
    gradMu = zeros(1,length(x));
    return
end

if nargin == 3 && isstruct(varargin{1})
    if ~isfield(varargin{1},'epsilon')&&length(fieldnames(varargin{1})) ~= 4 ...
            || isfield(varargin{1},'epsilon')&&length(fieldnames(varargin{1})) ~= 5
        error('Please input all optional parameters as structure')
    end
    theta = varargin{1}.theta;
    sigma = varargin{1}.sigma;
    ShiftType = varargin{1}.shifttype;
    NormWeighting = varargin{1}.NormWeighting;
else
    if nargin < 2
        error("Deflation requires the input of the deflated points, and the point to be evaluated")
    elseif nargin > 6
        error('Too many inputs')
    else
        if nargin > 2
            theta = varargin{1};
        else
            theta = 2; % default
        end
        if nargin > 3
            sigma = varargin{2};
        else
            sigma = 1; % default
        end
        if nargin > 4
            ShiftType = varargin{3};
        else
            ShiftType = "Multiple"; % default
        end
        if nargin > 5
            NormWeighting = varargin{4};
        else
            NormWeighting = speye(length(x)); % default
        end
    end
end

if length(x)~=size(Y,1)   
    % TODO: might be useful to allow length(x)<length(Y,1) and ignore any left over entries. 
    % Specifically of interest for "classic" IEP 
    error("The dimension of the current iterate and of deflated points must be the same.")
end


if strcmp(theta,"exp")  % strcmp needed because theta could be a double
    expi = exp(1 ./ sqrt(sum((NormWeighting*(x - Y)).^2,1)));

    switch ShiftType
        case "Single"
        MUs = expi -1;
        Mu = prod(MUs)+sigma;
        case "Multiple"
        MUs = expi + sigma - 1;
        Mu = prod(MUs);
        case "Scaled"
        MUs = (expi + sigma - 1).^(1/size(Y,2));
        Mu = prod(MUs);

    end
    if nargout>1
        gradMu = (-NormWeighting'*NormWeighting*sum(((Mu./MUs)./sum(abs(NormWeighting*(x-Y)).^2,1).^(3/2)).*expi.*(x-Y),2)).';
    end
else
    switch ShiftType
        case "Single"
            MUs = (1./sum(abs(NormWeighting*(x-Y)).^2,1).^(theta/2));
            Mu = sigma + prod(MUs);
        case "Multiple"
            MUs = sigma+ (1./sum(abs(NormWeighting*(x-Y)).^2,1).^(theta/2));
            Mu = prod(MUs);
        case "Scaled"
            MUs = (sigma+ (1./sum(abs(NormWeighting*(x-Y)).^2,1).^(theta/2))).^(1/size(Y,2));
            Mu = prod(MUs);
    end
        if nargout > 1
            gradMu = ((-theta*(NormWeighting'*NormWeighting))*sum(((Mu./MUs)./sum(abs(NormWeighting*(x-Y)).^2,1).^(1+theta/2)).*(x-Y),2)).';
        end
end
if ~exist('gradMu','var')
    gradMu = [];
end
if isinf(Mu)|| any(isnan(gradMu))
    if any(abs(Y - x)<1e-16)
        %  ME = MException('Deflation:AlreadyDeflated','The current iterate is a deflated point');
        % throwAsCaller(ME)
        warning('The current iterate is a deflated point')
    end
     % ME = MException('Deflation:NaNInf','Deflation operators calculated at the current point are Inf/NaN');
     %    throwAsCaller(ME)
        warning('Deflation operators calculated at the current point are Inf/NaN')
end
end

