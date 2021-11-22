function [k] = noyau(x, y)
    % linéaire
    %k = x'*y;

    % polynomial
    %k = (x'*y + 1)^2;

    % gaussien
    sigma = 5;
    k = exp(-norm(x-y)^2/(2*sigma^2));
end