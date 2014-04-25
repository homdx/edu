%%	Exercise 1 - kNN classifier

%%	Initialize
clear all; clf; clc;

%%	Part (a)   	Prepare the dataset and the classifier.
%	Generate data
n1 = 20; train_data_class1 = rand(n1,2);
n2 = 20; train_data_class2 = rand(n2,2) + ones(n2,2)*[1 0; 0 0];
train_data = [train_data_class1 ; train_data_class2];
train_label(1:n1,1) = 1;
train_label(n1+1:n1+n2,1) = 2;

%	Visualize the data
figure(1); clf; hold all; axis equal;
plot(train_data(1:n1,1), train_data(1:n1,2), 'r*');
plot(train_data(n1+1:n1+n2,1), train_data(n1+1:n1+n2,2), 'bo');
legend('Class 1', 'Class 2');
grid on;

%	Generate test points
n = 100;
test_data = [rand(n,2); rand(n,2) + ones(n,2)*[1 0; 0 0]];
test_label(1:n,1) = 1;
test_label(n+1:2*n,1) = 2;


%%	Part (b)   	Test the classifier.
%	Test the classifier.
k_values = [1, 3, 5, 7, 10, 15, 20];
errTrain = zeros(length(k_values), 1);
errTest = zeros(length(k_values), 1);
for i = 1:length(k_values)
	predTrain = knnClassify(train_data, train_label, train_data, k_values(i));
	predTest  = knnClassify(train_data, train_label, test_data,  k_values(i));
	errTrain(i) = loss01(predTrain, train_label);
	errTest(i)  = loss01(predTest, test_label);
end

%	Plot the training and the test errors.
figure(2); hold all;
plot(k_values, errTrain, 'r*:');
plot(k_values, errTest,  'b.-');

%	Plot the prediction for the best 'k'.
[~, bestIdx] = min(errTest);
bestk = k_values(bestIdx);
predTest  = knnClassify(train_data, train_label, test_data,  bestk);

figure(3); clf; hold all; axis equal;
pred_class1 = find(predTest == 1);
pred_class2 = find(predTest == 2);
plot(test_data(pred_class1, 1), test_data(pred_class1, 2), 'r*');
plot(test_data(pred_class2, 1), test_data(pred_class2, 2), 'bo');
plot([1 1], [0 1], 'k');
