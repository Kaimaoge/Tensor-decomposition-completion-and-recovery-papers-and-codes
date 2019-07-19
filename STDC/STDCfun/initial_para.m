function para = initial_para(kappa,omega,tau,gnns,pnns,maxitr,mode_dim,mode_PoM,mode_nse,VSet,Rate,Affinity,tsize)
%% Description
% VSet: the subsets encoding sub-manifolds, e.g., {[1,2],[3]} means that
%       the 1st and 2nd sub-manifolds are encoded simultaneously and the
%       3rd sub-manifold is encoded individually
% Rate: the downsampling rates of different Laplacian graphs
% Affinity: a cell array, where the i-th element denotes an affinity matrix
%           to encode the intra-factor relation w.r.t. the i-th tensor's
%           dimension

%% Initialization
para.kappa = kappa;               % the weight of MGE
para.omega = omega;               % (||Z||_F)^2 ~= omega*(ita*(||X-Zx1V1x2...xnVn||_F)^2)/2
para.tau = tau;                   % the 1st threshold for rank minimization
para.gnns = gnns;                 % # effective neighbors on manifold graphs
para.pnns = pnns;                 % # possible destinations on manifold graphs
para.maxitr = maxitr;             % maximal iteration of IALM algorithm
para.mode_dim = mode_dim;         % ignore the last dimension (true) or not (false)
para.mode_PoM = mode_PoM;         % update permutation matrices (true) or not (false)
para.mode_nse = mode_nse;         % ignore the observation noise (false) or not (true)

% convert graph indices into binary indicator 
num_g = numel(VSet);
N = numel(tsize);
if mode_dim, N=N-1; end;
for i = 1 : num_g
    id{i} = boolean(zeros(1,N));
    if sum(VSet{i}>N)>0
        error('Wrong index of manifold graph!');
    end
    id{i}(VSet{i}) = boolean(1);
end
para.VSet = id;
% ensure # downsampling rates == # graphs 
if numel(Rate)~=num_g
    disp('# graphs is incorrect...');
    disp('Correct it automatically...');
else
    Rate(num_g+1:end) = [];
    Rate(end+1:num_g) = 1;
end
para.Rate = Rate;
% confirm the affinity matrix
for i = 1 : N
    if isempty(Affinity{i}) || sum(abs(size(Affinity{i})-tsize(i)))~=0
        disp('Size of affinity matrix is incoorect...');
        disp('Correct it by default setting...');
        M = repmat(1:tsize(i),tsize(i),1);
        para.Affinity{i} = exp(-(M-M').^2/(2^2))-eye(tsize(i));
    else
        para.Affinity{i} = Affinity{i};
    end
end