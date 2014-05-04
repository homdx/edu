%%	Exercise 5 - Ridge regression


%%	Initialization
clear all; clc;
load('dataRidge.mat');


%%	Part (a)
figure(1);  clf;  hold all;
title('Ridge regression');  xlabel('x');  ylabel('y');

%	Plot the Training and Test points
scatter(x_train, y_train, 30, 'm', 'fill');	% Training points
scatter(x_test , y_test , 30, 'b', 'fill');	% Test points

%	Plot the fitted line
w = RidgeLLS(x_train, y_train, 0);
X = [min(x_train); max(x_train)];
y = [[1;1] X] * w;
line(X, y, 'LineStyle', '--', 'Color', 'r');

%	Predict and plot the test labels
y_test_pred = add1(x_test) * w;
scatter(x_test, y_test_pred, 40, 'r', 'fill');

%	Format
lgnd = char('y_{train}', 'y_{test}', 'y_{lin}', 'y_{test,lin}');
legend(lgnd);


%%	Part (b) - Linear least square with polynomial basis functions

%	Calculating the basis columns
powers = 1:15;

%	Non-linear (polynomial) features for x_train
X_train_poly = mapFeatures(x_train, powers);

w = RidgeLLS(X_train_poly, y_train, 0);

xx = -1.5:0.01:2.5;
Xx_poly = mapFeatures(xx, powers);
yy = add1(Xx_poly) * w;
plot(xx, yy, ':', 'Color', rand(1,3));
lgnd = char(lgnd, 'y_{poly}');
legend(lgnd);

%%
% All functions which may be represented or approximated by a 15 degree
% polynomial may be learnt by these basis functions, i.e.
% $\forall f \in P_{5}$


%%	Part (c) - Ridge regression

%	Plot the ridge regression function
lambda = [1e-4 1e-1 1e+1];
X_train_poly = mapFeatures(x_train, powers);
for i = 1:length(lambda)
	w = RidgeLLS(X_train_poly, y_train, lambda(i));
	y_train_ridge = add1(Xx_poly) * w;
	plot(xx, y_train_ridge, 'Color', rand(1,3));
	lgnd = char(lgnd, strcat('y_{ridge, \lambda=', num2str(lambda(i)), '}'));
end
legend(lgnd);

xlim([min([x_train; x_train]), max([x_train; x_train])]);
ylim([min([y_train; y_train]), max([y_train; y_train])]);

%	Prediction error
lambda = -15:1;
err = zeros(length(lambda), 1);
for i = 1:length(lambda)
	w = RidgeLLS(X_train_poly, y_train, 2 ^ lambda(i));
	err(i) = lossL2(y_train, add1(X_train_poly) * w);
end
[lambda',  2 .^ lambda',  err]
