function [X,info,itr,t] = HaLRTC(X,mark,alpha,ita,ita_rate,maxitr,N,X0)
tsize = size(X);
N0 = numel(tsize);
for i = 1 : N
    M{i} = zeros(tsize);
    Y{i} = zeros(tsize);
end
norm_x = norm(X(:));
tic
for itr = 1 : maxitr
    Xt = zeros(tsize);
    for i = 1 : N
        msize = circshift(tsize',N0+1-i);
        [Mi,D(i)] = Pro2TraceNorm(reshape(shiftdim(X+Y{i}/ita,i-1),tsize(i),[]),alpha(i)/ita);
        M{i} = shiftdim(reshape(Mi,msize'),N0+1-i);
        Xt = Xt+(M{i}-Y{i}/ita)/(N);
    end
    info.vrank(:,itr) = D;
    info.history(itr) = norm(X(mark)-Xt(mark),'fro');
    X(mark) = Xt(mark);
    
    for i = 1 : N
        Y{i} = Y{i}-ita*(M{i}-X);
        residual(i) = norm(X(:)-M{i}(:));
    end
    info.rse(itr) = norm(X(mark)-X0(mark),'fro')/norm(X0(:),'fro');
    info.residual(:,itr) = residual;
    
    if numel(find(residual>10^-4*norm_x))==0, break; end;
    ita = ita*ita_rate;
end
t=toc;
disp(['HaLRTC completed at ',int2str(itr),'-th iteration step within ',num2str(t),' seconds...']);