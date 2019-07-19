function [ B ] = generatef(V,U,j)
% U the cell matrix contain factor matrices of CP model
% B = [U[n-1],...,U[2],U[1],U[N],U[N-1],...U[n+1]]
dimT=length(U);
B=V;
 for i=1:dimT
     if i~=j
     B=nmodeproduct(B,U{i}',i);
     end
 end

end