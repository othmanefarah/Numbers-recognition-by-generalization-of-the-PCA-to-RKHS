function [k1, U] = PCA(Xc, PrecApprox)
[~,n] = size(Xc);
sigma = (Xc * Xc')/n;
[V,D] = eig(sigma);
[Val,ind] = sort(diag(D),'descend');
U = V(:,ind);
%C = U' * Xc;

y = 1 - sqrt(Val(2)/Val(1));
k1 = 2;
for i = 3 : n
    if (y < PrecApprox)
        y = 1 - sqrt(Val(i)/Val(1));
        k1 = i;
    else
        break;
    end
end

end