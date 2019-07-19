function [S, V, D] = MySVDtau(A, tau)
[m, n] = size(A);
if 10*m < n
    AAT = A*A';
    [S, V, D] = svd(AAT);
    V = diag(V) .^ 0.5;
    tol = max(size(A)) * eps(max(V));
    R = sum(V > max(tol, tau));

    %tol = min(size(A)) * eps(max(V));
    % R = sum(V > 1/tau)
        
%     tol = min(size(A)) * eps(max(V));
%     R = sum(V > tau)
    
    V = V(1:R);
    S = S(:,1:R);
    D = A'*S*diag(1./V);
    V = diag(V);
    return;
end
if m > 10*n
    [S, V, D] = MySVDtau(A', tau);
    mid = D;
    D = S;
    S = mid;
    return;
end
[S,V,D] = svd(A);
R = sum(diag(V) > tau);
S = S(:, 1:R);
V = V(1:R, 1:R);
D = D(:, 1:R);