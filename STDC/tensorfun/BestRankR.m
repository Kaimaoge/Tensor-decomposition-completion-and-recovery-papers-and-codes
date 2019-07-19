function [Ubasis Core Tensor k] = BestRankR(Tensor,maxitr,delta,Urank,mode,L,gmode)

% Initializaion
if mode==1
    Ubasis = HOSVD(Tensor,'N-1');
else
    Ubasis = HOSVD(Tensor,'N');
end
for d = 1 : ndims(Tensor)-mode
    Ubasis{1,d} = Ubasis{1,d}(:,1:Urank(d));
end

current_error = zeros(1,ndims(Tensor)-mode);
% Low rank approximation of best rank R
for k = 1 : maxitr
    last_error = current_error;
    for d1 = 1 : ndims(Tensor)-mode
        Z = Tensor;
        coef = 1;
        for d2 = 1 : ndims(Tensor)-mode
            if d1~=d2
                Z = TensorProduct(Z,Ubasis{1,d2}',d2);
                coef = coef*trace(Ubasis{1,d2}'*(L{d2}*L{d2})'*Ubasis{1,d2});
            end
        end
        if gmode==0 || gmode==1, coef = gmode; end;
        Z = shiftdim(Z,d1-1);
        Z = reshape(Z,size(Z,1),[]);
        [U S] = eig(Z*Z'-coef*(L{d1}*L{d1}'));
        [~,sidx] = sort(diag(S),'descend');
        U = U(:,sidx);
        current_error(d1) = norm(Ubasis{1,d1}'*U(:,1:Urank(d1)),'fro')/sqrt(length(Ubasis{1,d1}(:)));
        Ubasis{1,d1} = U(:,1:Urank(d1));
    end
    if length(find(current_error-(1-delta)*Urank<0)) == 0 | length(find(abs(current_error-last_error)>10^-7))==0
        break;
    end
end
clear Z;

% Rank-reduced approximation
for d = 1 : ndims(Tensor)-mode
    Tensor = TensorProduct(Tensor,Ubasis{1,d}',d);
end
Core = Tensor;
for d = 1 : ndims(Tensor)-mode
    Tensor = TensorProduct(Tensor,Ubasis{1,d},d);
end