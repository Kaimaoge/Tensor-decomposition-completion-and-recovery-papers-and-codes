function [Xm,A,R]=fused_tensor_completion_bywu(X,mark,varargin)
if ~exist('X','var'),    error('please input the sampled tensor.\n');    end
if ~exist('mark','var'),    error('please input the missing data indication mark.\n'); end
%%         default parameters   
           tsize=size(X);
           X(mark)=0;s=norm(X(:),'fro');theta=10^(-2.5)*s;  % the stop citeria
           alpha=0.2*ones(length(tsize),1);                 % parameters for l1 norm regularizer
           delta=0.8*ones(length(tsize),1);                 % parameters for smooth term
           beta=10*ones(length(tsize),1);                   % parameters for l2 norm regularizer, in general alpha/t=(1-beta)/t
           zeta=10*ones(length(tsize),1); 
           maxiter=1000;                                    % the maximum iteration
           facup=0.008;                                     % the parameter for updation of tensor rank
           Rmin=1;                                          % the minimum estimated tensor rank
           Rmax=max(round(max(tsize)/10),2);                % the maximum estimated tensor rank
           Rchange=1;                                
           for j=1:length(tsize)
           L=zeros(tsize(j)-1,tsize(j));                    % 1 -1
           for i1=1:tsize(j)-1                              % 0 1
               for j1=1:tsize(j)                            %      .
                 if i1==j1                                  %        1 -1
                    L(i1,j1)=1;                             %           1
                 end
                 if (j1-i1)==-1
                   L(i1,j1)=-1;
                 end
               end
           end
           Lp{j}=L'*L;  % the default smooth matrix (smooth variation=TV*TV') 
           Ls{j}=L;
           end
           type='normal';
%%         Read optional parameters
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for iv=1:2:(length(varargin)-1)
        switch upper(varargin{iv})
            case 'MAXITER',     maxiter=varargin{iv+1};
            case 'THETA',       theta=varargin{iv+1};
            case 'ALPHA',       alpha=varargin{iv+1};
            case 'BETA',        beta=varargin{iv+1};
            case 'DELTA',       delta=varargin{iv+1};
            case 'FACUP',       facup=varargin{iv+1};
            case 'RMIN',        Rmin=varargin{iv+1};
            case 'RMAX',        Rmax=varargin{iv+1};
            case 'LP',          Lp=varargin{iv+1};     % it should be noted that this input is a cell matrix
            case 'LS',          Ls=varargin{iv+1};     % it should be noted that this input is a cell matrix
            case 'ZETA',        zeta=varargin{iv+1};
            case 'TYPE',        type=varargin{iv+1};
            case 'RCHANGE',     Rchange=varargin{iv+1};
            otherwise
                error(['Unrecognized option: ',varargin{iv}]);
        end
    end
end

%% initiation
ma=max(X(:)); Xmo=ma/2+(ma/2)*rand(tsize);Xm=X; Xm(mark)=Xmo(mark); % initiation of missing entries
Xm(mark)=0;
R=Rmin;
for I=1:length(tsize)
 
A{I}=0.1*randn(tsize(I),R);  % We set the initiation factor matrices to be sufficient small

end

%% updation
for iter=1:maxiter
    factore=0;    % the updated matrix is initiated as 0
    tic
    for j=1:length(tsize)
        B=generate_cp_factor_matrx(A,j);
        V=unfolding(Xm,j);
        switch (type),
            case 'normal',
            [Anew]=fusedNeMF(V',B',A{j}',alpha(j),beta(j),delta(j),Lp{j}'); Anew=Anew';
            case 'nonnegative'
            [Anew]=fusedNeNMF(V',B',A{j}',alpha(j),beta(j),delta(j),Lp{j}',zeta(j),Ls{j}'); Anew=Anew';  
            case 'smooth'
            [Anew]=fusedNeMFs(V',B',A{j}',alpha(j),beta(j),delta(j),Lp{j}',zeta(j),Ls{j}'); Anew=Anew'; 
        end
            factore=factore+norm(Anew(:)-A{j}(:),'fro')/norm(A{j}(:),'fro');
        A{j}=Anew;    
    end
    Rold=R;
    if factore<=facup*length(tsize) && R<Rmax
          R=R+Rchange;
          for j=1:length(tsize)
             
            A{j}(:,Rold:R)=0.1*rand(tsize(j),R-Rold+1);
          end
    end
    
    %Xm1=double(ktensor(ones(R,1),A)); use ktensor with tensor-toolbox
    Xm1=folding(A{1}*generate_cp_factor_matrx(A,1),1,tsize);
    Xm=Xm1;
    Xm(~mark)=X(~mark);
  
    e=Xm-Xm1; e(mark)=0;stope(iter)=norm(e(:),'fro');
    if stope(iter)<=theta
        break;
    end
  
end
end