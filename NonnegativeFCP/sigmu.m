function [H1]=sigmu(H,alpha)

%H1= sign(H).* min( H/alpha, 1);
%hs=sign(H);
%H1(H1==inf)=hs(H1==inf);
H1 = H/alpha;
H1(H1>1) = 1;
H1(H1<-1) = -1;
end