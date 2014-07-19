%%	Exercise 3 - Generating samples from a Gaussian distribution


%%	Initialize
clear all; close all; clc;


%%	Part (b)
d = 3; n = 2000;
mu = [0, 0, 0];  Sigma = [2 0 0; 0 1 0; 0 0 4];
S1 = chol(Sigma);  X = repmat(mu,n,1) + randn(n,d)*S1;

figure(1);
plot3(X(:,1), X(:,2), X(:,3), 'r.');
grid on;


%%	Part (c)
[V, D] = eig(Sigma)


%%	Part (d)
A = [3 0; 0 1];
V = (1/sqrt(2)) * [1 1; -1 1];
d = 2; n = 400; mu = [0 0];
Sigma = V*A*V';
S1 = chol(Sigma); X = repmat(mu,n,1) + randn(n,d)*S1;

figure(2);  hold all;
line([-4 4], [4 -4]);
plot(X(:,1), X(:,2), 'r.');
grid on;
