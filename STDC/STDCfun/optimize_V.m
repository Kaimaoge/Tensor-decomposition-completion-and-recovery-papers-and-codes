function [V,P,rank_vi,pstr] = optimize_V(X,Y,Z,V,P,ita,tsize,vsize,para_ST)
%% Initialization
num_V = numel(V);
num_L = size(para_ST.H,1);
rank_vi = zeros(1,num_V);
pstr = '(';

%% Update by linearization
for i = 1 : num_V
    list = 1:num_V;
    list(i) = [];
    Bi = reshape(shiftdim(TensorChainProductT(Z,V,list),i-1),vsize(i),[]);
    Ci = -Bi*(reshape(shiftdim(Y,i-1),tsize(i),[])+ita*reshape(shiftdim(X,i-1),tsize(i),[]))';
    Bi = (Bi*Bi')*ita/2;
    
    Ai = zeros(tsize(i));
    for j = 1  : num_L
        if para_ST.VSet{j}(i)
            list = 1:num_V;
            list = list(para_ST.VSet{j});
            for k = 1 : numel(list)
                U{k} = V{list(k)}(:,P{list(k)})*para_ST.Ds{j,list(k)};
            end
            idx = find(list==i);
            list = 1:numel(list);
            list(idx) = [];
            if para_ST.Rate(j)==1, lidx=1; else lidx=idx; end;
            T = reshape(shiftdim(TensorChainProduct(para_ST.H{j,lidx},U,list),idx-1),tsize(i),[]);
            T = T*T';
            Ai = Ai+(T+T')/2;
        end
    end
    pidx = 1:tsize(i);
    pidx(P{i}) = pidx;
    At = Ai(pidx,pidx);
    At = At+At';
    Bt = Bi+Bi';
    sig = norm(At)+norm(Bt);
    grad = V{i}*At+Bt*V{i}+Ci;
    [Vi,r] = optimize_LRM(V{i},grad,sig,para_ST.alpha(i));
    
    V{i} = Vi;
    if para_ST.mode_PoM
        [pidx,info,cpu_time] = PermutationOnManifolds(V{i},At,tsize(i),para_ST.pnns,'euclidean',false);
        P{i} = pidx(P{i});
    end
    rank_vi(i) = r;
    pstr = [pstr,int2str(sum((1:tsize(i))'~=P{i})),','];
end
pstr(end) = ')';

function [optimum,r] = optimize_LRM(Mk,grad,sig,alpha)
%% M = arg min { alpha*|| M ||_* + <grad(Mk),M-Mk> + (sig/2)*|| M-Mk ||_F }
try
    [U,S,V] = svd(Mk-grad/(sig),'econ');
    S = diag(S);
    [S,r] = shrinkage(S,alpha/(sig));
    if r==0
        optimum = Mk;
    else
        optimum = U*diag(S)*V';
    end
catch
    disp('SVD fails to coverge...');
    optimum = Mk;
    r = 0;
end