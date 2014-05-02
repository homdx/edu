function [ X_poly ] = mapFeatures( x, powers )
%mapFeatures Map features to other powers.

X_poly = zeros(length(x), length(powers));
for i = powers
	X_poly(:,i) = x .^ i;
end

end
