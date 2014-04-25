%%	Exercise 4 - Decision boundary

%%	Initialization
clear all; clc;
[X, Y] = mixGaussian2d(100,0.4,0.6);

mu = zeros(2);
sigma = zeros(2);

mu(1,:) = mean(X(Y == 1));  sigma(1,:) = var(X(Y == 1));
mu(2,:) = mean(X(Y == 2));  sigma(2,:) = var(X(Y == 2));
