function Z = optimize_Z(V,X,Y,ita,gamma)
%% Closed form solution
% z = (ita*V*V^T + 2*gamma*I)^(-1) * V * (ita*x+y), where V=kron(V_1,...,V_n)
% -- 
% -- simplification
%
% let V_i*(V_i)^T = P_i*S_i*(P_i)^T, by spectral decomposition
% (ita*V*V^T + 2*gamma*I)^(-1) = 
%      kron(P_1,...,P_n) * (ita*kron(S_1,..S_n) + 2*gamma*I)^(-1) * kron((P_1)^T,...,(P_n)^T)
%     , where
%      (ita*kron(S_1,..S_n) + 2*gamma*I)^(-1) can be caculated efficiently

%%
wt = 1;
N = numel(V);

for i = ndims(X) : -1 : N+1
    wt = ones(size(X,i),1)*wt';
    wt = wt(:);
end
for i = N : -1 : 1
    [U,S] = eig(V{i}*V{i}');
    Q{i} = S;
    P{i} = U;
    wt = diag(S)*wt';
    wt = wt(:);
end
wt = reshape(ita*wt+2*gamma,size(X));  % diagonal vector of (ita*kron(S_1,..S_n) + 2*gamma*I)^(-1)
Z = TensorChainProduct(ita*X+Y,V,1:N);
Z = TensorChainProductT(Z,P,1:N);
Z = Z./wt;
Z = TensorChainProduct(Z,P,1:N);

