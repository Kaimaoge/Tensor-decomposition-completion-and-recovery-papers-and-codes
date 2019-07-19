function [Accuracy] = TensorRecognition_Rank1(U,Core,Data,Label,mode)
% Data: NxP
% Core: p1x...xpkxP
% mode: specify the mode of parameters for feature extraction

csize = size(Core);
Core = (reshape(Core,prod(csize)/csize(end),csize(end)))';
Coef0 = zeros(csize(mode),size(Data,1));
Coef0_Opt = zeros(csize(mode),size(Data,1));

CoefT = pinv(Core)*(Data');

for i = 1 : size(Data,1)
    coef0 = reshape(CoefT(:,i),csize(1:end-1));
    coef0 = reshape(coef0,csize(1:end-1));
    
    cvec0 = BestRank1(coef0,'',100,10^-5);
    cvec0_opt = BestRank1(coef0,'opt',100,10^-5);
    Coef0(:,i) = cvec0{mode};
    Coef0_Opt(:,i) = cvec0_opt{mode};
end

[ac0 ac0_per] = NN_Classifier(U{1}',Coef0,[1:size(U{1},1)]',Label,'cosine-t');
[ac0_opt ac0_opt_per] = NN_Classifier(U{1}',Coef0_Opt,[1:size(U{1},1)]',Label,'cosine-t');

Coef{1} = Coef0;
Coef{2} = Coef0_Opt;
Accuracy{1}.all = ac0;
Accuracy{2}.all = ac0_opt;
Accuracy{1}.per = ac0_per;
Accuracy{2}.per = ac0_opt_per;
