function [P,info,cpu_time] = PermutationOnManifolds(X,L,data_num,target_num,dType,sType)
%% initialization
% dimension of data/variables
target_num = min(target_num,data_num);
para_num = data_num*target_num;
constraint_num = 2*data_num;
% unknown variables
max_itr = 100;
P = (1:data_num)';
% constraint set
constraint_idx = zeros(2*para_num,1);
constraint_idx(1:para_num) = reshape(repmat((1:data_num)',1,target_num),1,[])';
para_idx = zeros(2*para_num,1);
para_idx(1:para_num) = 1:para_num;
constraint_result = ones(constraint_num,1);
% correlation/distance matrix
M_corr = X'*X;
M_dist = squareform(pdist(X',dType));
M_sqr2 = repmat(sum(X.^2,1),data_num,1);
M_sqr2 = M_sqr2+M_sqr2';
% optimization setting
opts = optimoptions('linprog','Algorithm','simplex','Display','off');
% algorithm information
cost = trace(X(:,P)*L*X(:,P)');
order = P;

L_diag = repmat(diag(L),1,data_num);
L_diff = L_diag+L_diag'-L-L';
%% main algorithm
% tic
for itr = 1 : max_itr
    % update the pairwise scores of every permutation
    if sType
    M_score = M_corr'*L+M_corr*L';
    M_score = M_score-repmat(diag(M_score),1,data_num);
    M_score = M_score+repmat(diag(M_corr),1,data_num).*L_diff+(M_corr+M_corr').*(L-L_diag)';
    else
    M_score = L*M_corr;
    M_score = (M_score+M_score')-(repmat(diag(M_score),1,data_num)+repmat(diag(M_score)',data_num,1));
    M_score = 2*M_score+(L_diff.*(M_sqr2-2*M_corr));
    end
    % find the target index
    [~,target_idx] = sort(M_dist-eye(data_num),2,'ascend');
    target_idx = target_idx(:,1:target_num);
    
    % define the permutation scores of source-to-target
    M_score = M_score(data_num*(target_idx-1)+repmat((1:data_num)',1,target_num));
    
    % define the constraint set for permutation
    tag = para_num;
    for i = 1 : data_num
        idx = find(target_idx==i);
        constraint_idx(tag+1:tag+numel(idx)) = i+data_num;
        para_idx(tag+1:tag+numel(idx)) = idx;
        tag = tag+numel(idx);
        constraint_result(i+data_num) = (numel(idx)~=0);
    end
    M_constraint = sparse(constraint_idx,para_idx,ones(2*para_num,1),constraint_num,para_num);
    
    % solve and update the permutation
    pidx = linprog(M_score(:),[],[],M_constraint,constraint_result,zeros(para_num,1),ones(para_num,1),[],opts);
    target_idx = target_idx';
    pidx = target_idx(reshape(pidx,data_num,[])'==1);
    P = P(pidx);% = P;
    
    M_corr = M_corr(pidx,pidx);
    M_sqr2 = M_sqr2(pidx,pidx);
    M_dist = M_dist(pidx,pidx);

    sidx = itr;
    % stopping criterion
    if sum(sum(abs(repmat(P,1,itr)-order),1)==0)>0
        [~,sidx] = min(cost); sidx = sidx(1); P = order(:,sidx);
        break; 
    end;
    cost = [cost;trace(L*M_corr)];
    order = [order,P];
%     disp([int2str(itr),'-th iteration of Manifold permutation (',int2str(sum((1:data_num)'~=P)),' data permuted, ',num2str(toc),' seconds)...']);
end
cpu_time = toc;
% disp([int2str(itr),'-th iteration of Manifold permutation (',int2str(sum((1:data_num)'~=P)),' data permuted, ',num2str(cpu_time),' seconds)...']);
info.cost = cost;
info.order = order;
info.sidx = sidx;
