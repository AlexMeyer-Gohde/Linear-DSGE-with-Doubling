function P=find_guess(A,B,C,maxeig)
    n=size(A,2);
    P=zeros(n,1);
    for j=1:n
        P(j)=find_pj(A(:,j),B(:,j),C(:,j),-maxeig,maxeig);
    end
    P=diag(P);
end

function pj=find_pj(aj,bj,cj,lb,ub)
    
    % Get coefficients:
    a = aj'*aj;
    b = 2*aj'*bj;
    c = bj'*bj + 2*aj'*cj;
    d = 2*bj'*cj;
    %e = cj'cj;
    
    % Solve cubic polynomial (only real solutions):
    coef = [4*a, 3*b, 2*c, d];
    pj=roots(coef);
    pj=real(pj(~imag(pj)));


    % Get minimizing pj in [lb,ub]
    pj=[lb;ub;pj((pj>=lb) & (pj<=ub))];
    [~,imin] = min(a*pj.^4 + b*pj.^3 + c*pj.^2 + d*pj);
    pj=pj(imin);

end
