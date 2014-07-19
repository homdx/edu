%%	Exercise 4 - PCA


%%	Initialize
clear all;  close all;  clc;


%%	Part (a)
%	See PCA.m


%%	Part (b)

%	Generate data
n = 500;  d = 2;
mu = [1 1];  Sigma = [2 -1; -1 2];
S1 = chol(Sigma); X = repmat(mu, n, 1) + randn(n, d) * S1;

%	Comparison of eigenvectors
[eigvecSigma, ~] = eig(Sigma)
pcaDir = PCA(X)
