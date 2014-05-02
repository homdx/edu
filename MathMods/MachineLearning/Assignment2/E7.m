%%	Exercise 7 - Least squares regression

%%	Initialize
clear all; clc;

sigma = 0.1;
[X_train, Y_train] = genLinData(50, sigma);
[X_test,  Y_test ] = genLinData(30, sigma);

%%	Part (a)

%	Description of the data
%	X contains uniformly distributed random numbers, s.t. X ? [-2, 1]
%	Y = 1 - 0.7*X + noise, where the noise varies on the variance ?.
figure(1);  clf;  hold all;  grid on;
scatter(X_train, Y_train, 'b*');
scatter(X_test, Y_test, 'g*');
xlabel('X');  ylabel('Y');
legend('Y_{train}', 'Y_{test}');

%	Preprocess - not required.

%%	Part (b)

%	Ploting the fitted line.
w = LLS(X_train, Y_train);
X = [min(X_train); max(X_train)];
y = [[1;1] X] * w;
line(X, y);
legend('Y_{train}', 'Y_{test}', 'Y_{pred}');

%	Predict and plot the test labels
Y_test_pred = [ones(size(X_test, 1), 1) X_test] * w;
scatter(X_test, Y_test_pred, 'r.');
legend('Y_{train}', 'Y_{test}', 'Y_{pred}', 'Y_{test,pred}');

%%	Part (c)

%	Average test error.
sigma = [0.01 0.1 0.4 0.9 1];
err = zeros(1, length(sigma));
for i = 1:length(sigma)
	for j = 1:10
		[X_train,Y_train] = genLinData(50, sigma(i));
		[X_test, Y_test]  = genLinData(30, sigma(i));
		w = LLS(X_train, Y_train);
		err(i) = err(i) + lossL2(Y_test, [ones(size(X_test, 1), 1) X_test] * w);
	end
end
err = err / 10;

figure(2); clf;	grid on;
plot(sigma, err);
xlabel('\sigma');  ylabel('Average error');
title('Error as a function of \sigma');

%	Sensitivity to outliers.
figure(1);
X_train = [X_train; 10];
Y_train = [Y_train; 10];
w = LLS(X_train, Y_train);
X = [min(X_train); max(X_train)];
y = [[1;1] X] * w;
line(X, y);
axis([-2 1, 0 2.5]);
legend('Y_{train}', 'Y_{test}', 'Y_{pred}', 'Y_{test,pred}', 'Y_{pred_{new}}');
%	The above plot shows that the model is highly sensitive to outliers.
