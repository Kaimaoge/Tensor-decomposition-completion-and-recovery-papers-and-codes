function [H]=fusedNeMFs(V,W0,H0,alpha,beta,delta,Lp,zeta,Ls)
if ~exist('V','var'),    error('please input the sample matrix.\n');    end
if ~exist('W0','var'),    error('please input the factor matrices W.\n'); end
if ~exist('H0','var'),    error('please input the factor matrices H.\n'); end
if ~exist('alpha','var'),    error('please input the sparse parameters alpha.\n'); end
if ~exist('beta','var'),    error('please input the smooth parameters beta.\n'); end
if ~exist('delta','var'),    error('please input the L2 norm parameters delta.\n'); end
if ~exist('Lp','var'),    error('please input the smooth/manifold matrix.\n'); end
if ~exist('zeta','var'),    error('please input the L1 smooth parameters delta.\n'); end
if ~exist('Ls','var'),    error('please input the smooth matrix.\n'); end
[m,n]=size(V);

tol=1e-5;
mu=1;
 
STOP_RULE = 3; 
W=W0; H=H0;
HVt=H*V'; HHt=H*H';
WtV=W'*V; WtW=W'*W;


%HIS.objf=sum(sum(WtW.*HHt))-2*sum(sum(WtV.*H))+alpha*sum(sum(abs(H)))+beta*sum(sum((L*H).^2));
LpC=norm(Lp); LsC=norm(Ls);
if ~issparse(WtW),
    L=norm(WtW+delta*eye(size(WtW)))+beta*LpC+alpha/mu+zeta/mu*LsC^2;	% Lipschitz constant
else
    L=norm(full(WtW+delta*eye(size(WtW))))+beta*LpC+alpha/mu+zeta/mu*LsC^2;
end
Grad=WtW*H-WtV+alpha/mu*sigmu(H,alpha/mu)+beta*H*Lp+delta*H+zeta/mu*sigmu(H*Ls,mu)*Ls'; 

alpha1=1;
Z=H0;
ITER_MAX=5;      % maximum inner iteration number (Default)
ITER_MIN=2;
for iter=1:ITER_MAX,
    H0=H;
    H=Z-Grad/(L+1e-10);    % Calculate sequence 'Y'
    alpha2=0.5*(1+sqrt(1+4*alpha1^2));
    Z=H+((alpha1-1)/alpha2)*(H-H0);
    alpha1=alpha2;
    Grad=WtW*Z-WtV+alpha/mu*sigmu(Z,mu)+beta*Z*Lp+delta*Z+zeta/mu*sigmu(Z*Ls,mu)*Ls'; 
    % Stopping criteria
    if iter>=ITER_MIN,
        pgn=GetStopCriterion(STOP_RULE,Z,Grad);
        if pgn<=tol,
            break;
        end
    end
end



end