function [ B ] = generate_cp_factor_matrx(U,j)
% U the cell matrix contain factor matrices of CP model
% B = [U[n-1],...,U[2],U[1],U[N],U[N-1],...U[n+1]]
dimT=length(U);
if j==1 || j==dimT
        for k=1:dimT
            Uk{k}=U{dimT+1-k};
        end
        Uk(dimT+1-j)=[];
else
     for k=1:j-1
         Uk{k}=U{j-k};
     end
     for k=j:dimT-1
         Uk{k}=U{dimT+j-k};
     end    
end
    B=(kr(Uk))';
end