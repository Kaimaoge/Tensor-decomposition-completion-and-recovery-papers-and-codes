clc
clear all
for ratio=1:19
    
S=rand(15,15,15);
X=rand(100,15);
Y=rand(100,15);
Z=rand(100,15);

SX=nmodeproduct(S,X,1);
SXY=nmodeproduct(SX,Y,2);
L=nmodeproduct(SXY,Z,3);
max_w=max(reshape(L,100*100*100,1));
L=L/max_w;
S=zeros(100,100,100);
for m=1:100
    for n=1:100
        for k=1:100
            S(m,n,k)=0.2*(rand(1))-0.1;
        end
    end
end
S(randperm(100*100*100)<=100*100*100*(1-0.5))=0;
S1=zeros(100,100,100);
for m=1:100
    for n=1:100
        for k=1:100
            S1(m,n,k)=2*(rand(1))-1;
        end
    end
end
S1(randperm(100*100*100)<=100*100*100*(1-0.3))=0;
S2=zeros(100,100,100);
for m=1:100
    for n=1:100
        for k=1:100
            S2(m,n,k)=20*(rand(1))-10;
        end
    end
end
S2(randperm(100*100*100)<=100*100*100*1-0.1)=0;

W=ones(100,100,100);
W(randperm(100*100*100)<=(100*100*100*0.05*ratio))=0;



M=L+S+S1+S2;
MW=M.*W;
Omega=logical(W);
alpha=[1,1,1e-3];
maxIter = 1000;
epsilon = 1e-10;
factor = 2;
C =  0.6;
L0 = 1e-5;
[Take, errList_H] = FaLRTCnr(...
     MW,...                    % a tensor whose elements in Omega are used for estimating missing value
     Omega,...          % the index set indicating the obeserved elements
     alpha,...             % the coefficient of the objective function,  i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_* 
     factor,...             % control the decay speed of the relaxation parameter \mu, i.e., increase \mu by \mu = O(1 / k^factor) where k denotes the kth iteration. A reasonable "factor" could be in the range (0, 4].  
     L0,...                   % the initial step size parameter, i.e., stepsize = 1 / L;
     C,...                     % the shrinkage ratio in the range (0.5, 1). When "L" does not pass the descent enough test, we update L = L / C; 
     maxIter,...        % the maximum iterations
     epsilon...          % the tolerance of the relative difference of outputs of two neighbor iterations 
    );


W1=double(1-W);
ee=(L-Take).*W1;
error(ratio)=norm(unfolding(ee,1),'fro')/norm(unfolding(L.*W1,1),'fro');
end
error=error';