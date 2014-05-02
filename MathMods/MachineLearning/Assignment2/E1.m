%%	Exercise 1 - Text classification

%%	Initialization
clear all; clc;
load('20Newsgroup.mat');
% n = 40;

%%	Training data set
trIdx6 = find(y_train == 6);  rand_train_6 = randperm(length(trIdx6));
trIdx8 = find(y_train == 8);  rand_train_8 = randperm(length(trIdx8));
% trList = [trIdx6(rand_train_6(1:n));
%           trIdx8(rand_train_8(1:n))];
trList = [trIdx6;
          trIdx8];
x_train_6_8 = x_train(trList,:);  y_train_6_8 = y_train(trList);

%%	Test data set
x_test = x_test(:, 1:size(x_train, 2));    % Clean the data
tsIdx6 = find(y_test == 6);  rand_test_6 = randperm(length(tsIdx6));
tsIdx8 = find(y_test == 8);  rand_test_8 = randperm(length(tsIdx8));
% tsList = [tsIdx6(rand_test_6(1:n));
%           tsIdx8(rand_test_8(1:n))];
tsList = [tsIdx6;
          tsIdx8];
x_test_6_8 = x_test(tsList,:);  y_test_6_8 = y_test(tsList);

%%	Use kNN Classifier
k = 1:2:9;
err_train = zeros(length(k), 1);
err_test = zeros(length(k), 1);
for i = 1:length(k)
	pred_train = knnClassify(x_train_6_8, y_train_6_8, x_train_6_8, k(i));
	pred_test  = knnClassify(x_train_6_8, y_train_6_8, x_test_6_8,  k(i));
	err_train(i) = loss01(pred_train, y_train_6_8);
	err_test(i)  = loss01(pred_test,  y_test_6_8);
end

%%	Plot the training and the test errors.
figure(1); hold all; grid on;
plot(k, err_train, 'b*-');
plot(k, err_test,  'ro-');
legend('Training error', 'Test error');
