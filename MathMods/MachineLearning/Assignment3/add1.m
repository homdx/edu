function [ X ] = add1( X )
%add1 Adds a column of 1s in front of the matrix

X = [ones(size(X,1), 1) X];

end
