function [ err ] = loss01( y_pred, y )
%loss01 Gives an estimate of the error in prediction
%   Uses the L1 norm.

err = mean(abs(sign(y_pred - y)));

end
