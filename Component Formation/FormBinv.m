function [Binv,B,CholR]= FormBinv(A)
%FORMBINV Forms inverse of B.
%   Binv = FORMBINV(A) calculates the inverse of the matrix B, as used in the
%   Lift and Projection method, where A is a cell array of the basis matrices
%   of A(x).
%
%   [Binv,B]= FormBinv(A) also returns the matrix B.
%
%   [Binv,B,CholR]= FormBinv(A) also returns the Choelesky factor of B.
l = length(A);
B = zeros(l);
for i = 1:(l)
    B(i,i) = sum(sum(A{i}'.*A{i}));
    for j = i+1:(length(A))
        B(i,j) = A{j}(:)'*A{i}(:);
        B(j,i) = B(i,j);
    end
end
B = sparse(B);
Binv = inv(B);
if nargout>2
    CholR = chol(B);
end
end