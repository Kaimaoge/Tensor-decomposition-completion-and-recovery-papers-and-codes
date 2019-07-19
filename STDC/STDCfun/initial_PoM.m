function [P,info] = initial_PoM(X,nn,num,N)
tsize = size(X);
for i = 1 : N
    T = reshape(shiftdim(X,i-1),tsize(i),[])';
    [~,S,V] = svd(T,'econ');
    T = S*V';
    W = repmat(1:tsize(i),tsize(i),1);
    W = exp(-(abs(W-W')/nn).^2);
    L = diag(sum(W,2))-W;
    [P{i},info{i},cpu_time] = PermutationOnManifolds(T,L,tsize(i),num,'euclidean',false);
end