function H_lambda = FormH_Lambda(Q,A,D,ev)
    %Form S
    m = length(ev); l=length(A);    
    H_lambda = zeros(l,l,m);
    QAQ = cell(1,l);
    for i = 1:l
        QAQ{i} = Q'*A{i}*Q;
        QAQ{i} =QAQ{i}(1:m,1:m);
    end
    DD=D'-D;
    DD(abs(DD)<1e-15) = Inf;
        for j=1:l
            for k = 1:l
                H_lambda(k,j,:) = (2*sum(QAQ{k}.*QAQ{j}./DD));
            end
        end
end
 

