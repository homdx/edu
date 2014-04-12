%%	Group 1 - Assignment 1
%	Anirudh Gupta		-	6607652
%	Maria Panoukidou	-	6607946
%	Sudip Sinha		-	6609779

%%	Initialize
clear all; clf; clc;
digits = [2 3];	% Change this to [3 8] for the last part of the exercise.


%%	Part (a)   	Prepare the training data.

%	Train
load('usps_train.mat');
trList = find(train_label == digits(1) | train_label == digits(2));
x_train = double(train_data(trList,:));
y_train = double(train_label(trList));

%	Test
load('usps_test.mat');
tsList = find(test_label == digits(1) | test_label == digits(2));
x_test = double(test_data(tsList,:));
y_test = double(test_label(tsList));

%	Visualize the training example number 155.
% dig = reshape(train_data(155,:),16,16);
% imagesc(dig);
% colormap('gray');


%%	Part (b)   	Evaluate the performance of the classifier.

%	Test the classifier.
k_values = [1 3 5 7 10 15];
errTrain = zeros(length(k_values), 1);
errTest  = zeros(length(k_values), 1);
for i = 1:length(k_values)
	predTrain = knnClassify(x_train, y_train, x_train, k_values(i));
	predTest  = knnClassify(x_train, y_train, x_test,  k_values(i));
	errTrain(i) = loss01(predTrain, y_train);
	errTest(i)  = loss01(predTest,  y_test );
end

%	Plot the training and the test errors.
% figure(2);
hold all;
% grid on;
plot(k_values, errTrain, 'r*:');
plot(k_values, errTest,  'b.-');
legend('Train', 'Test');
