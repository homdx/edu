function [ w ] = RidgeLLS( X, y, lambda )
%RidgeLLS Linear Least Squares with L2 regularization
%   w = argmin((1/m)*(||y - Xw||_2)^2 + lambda*||w||_2)

X = add1(X);
w = (X' * X + lambda * size(X, 1) * eye(size(X, 2))) \ (X' * y);

end
