function [Ubasis Core Sigma] = HOSVD(Tensor,mode)

switch mode
    case 'N-1'
        for d = 1 : ndims(Tensor)-1
            T = shiftdim(Tensor,d-1);
            T = reshape(T,size(T,1),prod(size(T))/size(T,1));
            [Ubasis{d} S] = svd(T);
            Sigma{d} = diag(S);
        end
        Core = Tensor;
        for d = 1 : ndims(Tensor)-1
             Core = TensorProduct(Core,(Ubasis{d})',d);
        end
    case 'N'
        for d = 1 : ndims(Tensor)
            T = shiftdim(Tensor,d-1);
            T = reshape(T,size(T,1),prod(size(T))/size(T,1));
            [Ubasis{d} S] = svd(T);
            Sigma{d} = diag(S);
        end
        Core = Tensor;
        for d = 1 : ndims(Tensor)
             Core = TensorProduct(Core,(Ubasis{d})',d);
        end
end
