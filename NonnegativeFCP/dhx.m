clc
clear all
close all
%% understanding the spectral of traffic flow
%load('traffic_volume.mat')
%load('mantra.mat')

load('T.mat')
w1{1} = ones(21,21);
w1{2} = ones(21,21);
w1{3} = ones(24,24);

tsize = size(T);
om=ones(tsize);
r=0;
om(randperm(prod(tsize))<=r*prod(tsize))=0;
v_missing=T.*om;
[PP2,A,R]=fused_tensor_completion_bywu(T,~om,'RMIN',10,'BETA',[1,1,1],'RMAX',10,'ALPHA',[0.15,0.15,0.15],'DELTA',[0,0,0],'ZETA',[0,0,0],'TYPE','nonnegative','LP',w1,'RCHANGE',20,'FACUP',0.05);

figure 
subplot(3,1,1)
plot(A{1}(:,2))
subplot(3,1,2)
plot(A{2}(:,2))
subplot(3,1,3)
plot(A{3}(:,2))
