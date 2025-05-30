function H = FormH_Lambda(Q,A,D,evORm)
%FORMJLAMBDA Hessian matrix of Lambda
%
%FORMJLAMBDA(Q,A) is the Hessian matrix of Lambda with respect to x, where
%Q is the eigenvectors of the current Hamiltonian estimate, A is the 
% cell tructure of basis matrices, D is the eigenvalues of the current 
% Hamiltonian estimate, and evORm is either the vector if prescribed
% eigenvalues or the length of this vector.
  
if isscalar(evORm)
    m = evORm;
else
    m = length(evORm);
end
l=length(A);
H = zeros(l,l,m);
QAQ = cell(1,l);
for i = 1:l
    QAQ{i} = Q(:,1:m)'*A{i}*Q(:,1:m);
end
if ~isvector(D)
    D  =diag(D);
end
DD=D(1:m)'-D;
DD(abs(DD)<1e-15) = Inf;
for j=1:l
    for k = 1:l
        H(k,j,:) = real(2*sum(QAQ{k}.*QAQ{j}./DD));
    end
end