function [ w ] = LLS( X, y )
%LLS Linear Least Squares
%   w = argmin(||y - Xw||)

X = [ones(size(X, 1), 1) X];
w = pinv(X' * X) * (X' * y);

end
