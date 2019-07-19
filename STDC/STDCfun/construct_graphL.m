function H = construct_graphL(tsize,VSet,Rate,nbs,Affinity)
%% Initialization
% # sub-graphs
M = numel(VSet);
% range of neighboring indeces 
nbs = -nbs:nbs;
H = [];

%% Main algorithm
for m = 1 : M
    % parameters in the m-th graph
    
    % # sub-manifolds
    mnum = sum(VSet{m});
    % possible neighborhood
    NSet = nchoosek(repmat(nbs,1,mnum),mnum);
    NSet(sum(abs(NSet),2)==0,:) = [];
    % # neighbors
    Nnum = size(NSet,1);
    
    % construction of the (downsampled) m-th graph
    if Rate(m)==1, mk = 1; else mk = mnum; end;
    aidx = 1:numel(tsize);
    aidx = aidx(VSet{m});
    for k = 1 : mnum
        Am{k} = imresize(Affinity{aidx(k)},round(tsize(aidx(k))*Rate(m))*ones(1,2));
    end
    for k = 1 : mk
        ridx = ones(1,mnum)*Rate(m);
        ridx(k) = 1;
        As = Am;
        As{k} = Affinity{aidx(k)};
        msize = round(tsize(aidx).*ridx);
        N = prod(msize);
        % fill-in Laplacian matrix
        xidx = zeros(N*Nnum,1);
        yidx = xidx;
        zidx = xidx;
        len = 0;
        for j = 1 : N
            jidx = my_ind2sub(msize,j)';
            nidx = NSet+repmat(jidx,[Nnum,1]);
            nidx(sum(nidx<1,2)>0 | sum(nidx>repmat(msize,[Nnum,1]),2)>0,:) = [];
            nnum = size(nidx,1);
            edgeW = zeros([nnum,mnum]);
            for i = 1 : mnum
                edgeW(:,i) = As{i}((nidx(:,i)-1)*size(As{i},1)+jidx(i));
            end
            xidx(len+1:len+nnum) = j;
            yidx(len+1:len+nnum) = my_sub2ind(msize,nidx);
            zidx(len+1:len+nnum) = sum(edgeW,2);
            len = len+nnum;
        end
        L = sparse(xidx(1:len),yidx(1:len),zidx(1:len),N,N);
        L = (L+L')/2;
        % factorize Laplacian matrix
%         L = full(L);
        H{m,k} = full(cholcov(diag(sum(L,2))-L));
        H{m,k} = reshape(H{m,k}',[msize,size(H{m,k},1)]);
    end
    disp(['Construction of the ',int2str(m),'-th sub-graph completed...']);
end
