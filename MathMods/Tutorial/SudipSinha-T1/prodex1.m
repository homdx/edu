function prod = prodex1(x)
%prodex1 Product of a matrix or vector

prod = 1;
l = length(x);
for idx = 1:l
	prod = prod * x(idx);
end
end