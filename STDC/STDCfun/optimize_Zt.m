function [Z er_z gdnorm] = optimize_Zt(Z,V,X,Y,ita,gamma)
%% Initialization

maxitr = 50;
er_z = zeros(maxitr,1);
gdnorm = zeros(maxitr,1);

N = numel(V);
for i = 1 : N
    VVt{i} = V{i}*V{i}';
end
r = TensorChainProductT(Z,V,1:N)-(X+Y/ita);
r = -ita*TensorChainProduct(r,V,1:N) - 2*gamma*Z;
p = r;
rsold = sum(r(:).^2);
tol = 10^-3*sqrt(rsold);

%% Update by conjugate gradient method
for i = 1 : maxitr
    Ap = ita*TensorChainProduct(p,VVt,1:N) + 2*gamma*p;
    alpha = p.*Ap;
    alpha = rsold/sum(alpha(:));
    Z = Z+alpha*p;
    r = r-alpha*Ap;
    rsnew = sum(r(:).^2);
    
    er_z(i) = norm(alpha*p(:));
    gdnorm(i) = sqrt(rsnew);
    if gdnorm(i)<tol, break; end;
    
    p = r+(rsnew/rsold)*p;
    rsold = rsnew;
end
