%%	Group 1 - Assignment 1
%	Anirudh Gupta     -  6607652
%	Maria Panoukidou  -  6607946
%	Sudip Sinha       -  6609779

function [ err ] = loss01( y_pred, y )
%loss01 Gives an estimate of the error in prediction
%   Uses the L1 norm.

err = mean(y_pred ~= y);

end
