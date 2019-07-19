function [X r] = shrinkage(X,t)
idx1 = X>t;
idx2 = X<-t;
idx3 = X>=-t & X<=t;
X(idx1) = X(idx1)-t;
X(idx2) = X(idx2)+t;
X(idx3) = 0;
r = sum(X~=0);