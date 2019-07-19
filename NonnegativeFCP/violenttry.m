clc
clear all
close all
load('missing.mat')
load('true.mat')
for i=1:10
    for j=1:10
        for k=1:10
             beta=[10*k,10*k,5*k];
             alpha=[0.1*i,0.1*i,0.1*i];
             delta=1-alpha;
             zeta=[0.3*j,0.3*j,0.1*j];
[Xm,A,R]=fused_tensor_completion_bywu(X,mark,'TYPE','smooth','RMAX',60,'BETA',beta,'ALPHA',alpha,'DELTA',delta,'ZETA',zeta);
P=PSNR(Xm,X);
if P>24.1
    break;
end
        end
    end
end