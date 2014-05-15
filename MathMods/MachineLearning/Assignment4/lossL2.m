function [ err ] = lossL2( y, y_pred )
%lossL2 L2 loss function.

err = mean(norm(y - y_pred));

end

