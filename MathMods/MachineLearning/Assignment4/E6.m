%%	Exercise 6 - Cross validation in ridge regression


%%	Initialize
clear all; clc;
load('dataRidge.mat');
jj = -15:8;
lambda = 2 .^ (-jj);
meanErr = zeros(1, length(jj));


%%	Partition and compute errors
CVO = cvpartition(y_train,'kfold',10);
for i = 1:length(lambda)
    err = zeros(CVO.NumTestSets,1);
        for j = 1:CVO.NumTestSets
            trIdx = CVO.training(j);
            teIdx = CVO.test(j);
            w = RidgeLLS([x_train(trIdx,:) ones(CVO.TrainSize(j),1)], y_train(trIdx,:) ,lambda(i));
            Y_test_pred = polyval(w,x_train(teIdx,:));
            err(j) = err(j) + lossL2(y_train(teIdx,:),Y_test_pred);
        end
        meanErr(i) = sum(err)/CVO.NumTestSets;
end


%%	Display
plot(meanErr);
[~,mIdx] = min(meanErr);
fprintf('The best performance is for lambda = %f, corresponding to j = %i.\n', lambda(mIdx), jj(mIdx));
