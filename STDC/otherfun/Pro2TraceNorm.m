function [X, n] = Pro2TraceNorm(Z, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min: 1/2*||Z-X||^2 + ||X||_tr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[S, V, D] = MySVD(Z);
[S, V, D] = MySVDtau(Z, tau);
V = max(diag(V) - tau, 0);
n = sum(V > 0);
X = S(:, 1:n) * diag(V(1:n)) * D(:, 1:n)';




