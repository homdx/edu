function [ coeff ] = PCA( X )
%PCA    Principal Component Analysis

X_centered = X - repmat(mean(X), [size(X, 1) 1]);
Covariance = X_centered' * X_centered;
% [eigvec, eigval] = eig(Covariance);
% [~, order] = sort(diag(eigval), 'descend');
% coeff = eigvec(:, order);
[coeff, ~] = eig(Covariance);

end
