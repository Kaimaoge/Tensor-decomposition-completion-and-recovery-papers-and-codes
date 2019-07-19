function [V,Core,info,X,itr] = STDC(X,mark,para_ST,Xg,src_type,B_SVD)
%% Simultaneous Tensor Decomposition and Completion
% Input:
%  -- X: an input N order tensor object
%  -- mark: boolean index of missing entries (1 indicates missing; 0 otherwise)
%  -- para_ST: parameters used for the proposed STDC algorithm
%  -- mode: determining whether the N-th submanifold (V_N, usually ignored in multilinear model analysis) is computed (set as 0) or not (set as 1)
%  -- Xg: the ground truth of X; if not available, just using X instead

%% Initialization
tic;
tsize = size(X);
% initialize manifold graphs
if ~isfield(para_ST,'VSet')
    para_ST.H = [];
else
    para_ST.H = construct_graphL(tsize,para_ST.VSet,para_ST.Rate,para_ST.gnns,para_ST.Affinity);
    for i = 1 : size(para_ST.H,1)
        if size(para_ST.H,2)==1 || numel(para_ST.H{i,2})==0
            for j = 1 : numel(tsize)
                para_ST.Ds{i,j} = eye(tsize(j));
            end
        else
            for j = 1 : numel(tsize)
                A = randn(tsize(j)); 
                B = imresize(A,[tsize(j) round(tsize(j)*para_ST.Rate(i))],'bilinear');
                para_ST.Ds{i,j} = A\B;
            end
        end
    end
end
% initilize factor matrices & permutation matrices (w.r.t the modes of tensor dimension & PoM)
vsize = tsize;
N = numel(tsize);
if para_ST.mode_dim, N = N-1; end;
if para_ST.mode_PoM
    randn('seed',1);
    for i = 1 : N
        V{i} = randn(tsize(i));
        V{i} = V{i}/norm(V{i});
        vsize(i) = size(V{i},1);
    end
    Xt = HaLRTC(X,mark,ones(N,1),10^-2,1.1,100,N,Xg);
    P = initial_PoM(Xt,para_ST.gnns^2,para_ST.pnns,N);
else
    for i = 1 : N
        V{i} = eye(tsize(i));
        vsize(i) = size(V{i},1);
        P{i} = (1:tsize(i))';
    end
end
% initialize core tensor and augmented multiplier
Z = X;
Y = zeros(tsize);
% initialize algorithm parameters
norm_gt = norm(Xg(:));
norm_x = norm(X(:));
para_ST.alpha = ones(N,1);
para_ST.gamma = (para_ST.omega/para_ST.tau)/(norm_x^2);
xxt = reshape(X,tsize(1),[]);
xxt = norm(xxt*xxt');
ita = 1/(para_ST.tau*xxt);
ct = zeros(1,N);
for i = 1 : size(para_ST.H,1)
    ct = ct+double(para_ST.VSet{i});
end
for i = 1 : size(para_ST.H,1)
    list = 1 : N;
    list = list(para_ST.VSet{i});
    for j = 1 : size(para_ST.H,2)
        U{j} = para_ST.Ds{i,list(j)};
    end
    for j = 1 : size(para_ST.H,2)
        llt = reshape(TensorChainProduct(para_ST.H{i,j},U,[1:j-1 j+1:size(para_ST.H,2)]),tsize(list(j)),[]);
        llt = norm(llt*llt');
        para_ST.H{i,j} = para_ST.H{i,j}*para_ST.kappa*sqrt(ita*xxt/(2*llt*ct(list(j))));
    end
end
% message
disp(['Finish the initialization of all parameters within ',num2str(toc),' seconds...']);
disp('------------------------------------------------------------------------------');
disp('--                          Start STDC algorithm..                          --');
disp('------------------------------------------------------------------------------');
%% Main algorithm
lambda = ita*(1.1^(para_ST.maxitr))/2;
tic;
switch src_type
    case 'image'
        figure('Position',get(0,'ScreenSize'));
        subplot(1,3,1);imshow(X);title('incomplete tensor');
    case 'CMU'
        figure('Position',get(0,'ScreenSize'));
        subplot(1,3,1);imshow(reshape(permute(reshape(TensorProduct(X(1,:,:,:),B_SVD,4),[11,21,32,32]),[3,2,4,1]),32*21,[]));title('incomplete tensor (1st subject)');
    otherwise
end
%     aviname = input('input the file name for avi: ','s');
%     aviobj=avifile(aviname);   %定义一个视频文件用来存动画
%     aviobj.quality=60;
%     aviobj.Fps=5;
for itr = 1 : para_ST.maxitr
    % update V1,...,Vn
    [V,P,rank_vi,pstr] = optimize_V(X,Y,Z,V,P,ita,tsize,vsize,para_ST);
    % update Z
    Z = optimize_Z(V,X,Y,ita,para_ST.gamma);
    % update X
    Xt = TensorChainProductT(Z,V,1:numel(V));
    X(mark) = Xt(mark)-Y(mark)/ita;
    if para_ST.mode_nse
        X(~mark) = ((ita*Xt(~mark)-Y(~mark))+lambda*Xg(~mark))/(ita+lambda);
    end
    residual = (norm(X(:)-Xt(:)))/norm(Xt(:));
    % update Y
    Y = Y+ita*(X-Xt);
    % assessment
    info.rse(itr) = norm(X(mark)-Xg(mark))/norm_gt;
    info.rank_vi(:,itr) = rank_vi;
    info.residual(:,itr) = residual;
    % display
    disp_t = ['STDC completed at ',int2str(itr),'-th iteration step within ',num2str(toc),' seconds...(',pstr,' data are permuted)'];
    switch src_type
        case 'image'
            subplot(1,3,2);plot(info.rse);axis([0,para_ST.maxitr,0,inf]);title('# iterations vs. RSEs');
            subplot(1,3,3);imshow(X);title('completed tensor');
            axes('position',[0,0,1,1],'visible','off');
            text(0.05,0.98,['### ',disp_t]);
            pause(0.1);
        case 'CMU'
            subplot(1,3,2);plot(info.rse);axis([0,para_ST.maxitr,0,inf]);title('# iterations vs. RSEs');
            subplot(1,3,3);imshow(reshape(permute(reshape(TensorProduct(X(1,:,:,:),B_SVD,4),[11,21,32,32]),[3,2,4,1]),32*21,[]));title('completed tensor (1st subject)');
            axes('position',[0,0,1,1],'visible','off');
            text(0.05,0.98,['### ',disp_t]);
            pause(0.1);
           
        otherwise
    end
%     disp(disp_t);
%             frame=getframe(gca);   %把图像存入视频文件中
%             im=frame2im(frame);
%             aviobj=addframe(aviobj,im);
    % stopping criterion
%    if residual<1e-2, break; end;
    ita = ita*1.1;
end
% aviobj=close(aviobj);
Core = Z;
for i = 1 : N
    V{i} = V{i}';
end
info.P = P;
