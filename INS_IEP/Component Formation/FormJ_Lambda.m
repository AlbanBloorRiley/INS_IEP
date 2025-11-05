function J = FormJ_Lambda(Q,A)
%FORMJLAMBDA Jacobian matrix of Lambda
%
%FORMJLAMBDA(Q,A) is the Jacobian matrix of Lambda with respect to x, where
%Q is the eigenvectors of the current Hamiltonian estimate, and A is the 
% cell tructure of basis matrices.

l = length(A);
J = zeros(size(Q,2),l);
for k = 1:l
    J(:,k) =real(sum((Q.'*A{k}).*Q',2));
end

