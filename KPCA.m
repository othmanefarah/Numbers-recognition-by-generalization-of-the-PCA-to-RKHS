function [k2, K, alpha, Y, V0] = KPCA(Xc, PrecApprox)
[~,n] = size(Xc);
K = zeros(n, n);
for i = 1 : n
    for j = 1 : n
        K(i,j) = noyau(Xc(:,i), Xc(:,j)); 
    end
end
Kc = K - ones(n,n)*K/n - K*ones(n,n)/n + ones(n,n)*K*ones(n,n)/(n^2);
[V, D] = eig(Kc);
[Val, ind] = sort(diag(D), 'descend');
U = V(:, ind);

y = 1 - sqrt(Val(2)/Val(1));
k2 = 2;
for i = 3 : n
    if (y < PrecApprox)
        y = 1 - sqrt(Val(i)/Val(1));
        k2 = i;
    else
        break;
    end
end

alpha = zeros(n, k2);
for i = 1 : k2
    alpha(:,i) = U(:,i)/sqrt(Val(i));
end

Y = alpha' * K';
V0 = U(:, 1:k2);
end