function B = Form_Gram_B(A)
%Forms the Gram matrix B used in LP method
%
%Only requires the basis matrices - A - as input
if iscell(A)
    l = length(A);
    B = zeros(l);
    for i = 1:l
        for j = 1:l
            B(i,j) = sum(sum(A{j}'.*A{i}));
        end
    end
else
    l = size(A,3);
    B = zeros(l);
    for i = 1:l
        for j = 1:l
            B(i,j) = sum(sum(A(:,:,j)'.*A(:,:,i)));
        end
    end
end
end
