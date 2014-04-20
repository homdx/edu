function sum = sumex1(x)
%sumex1 Summing up a matrix or vector

sum = 0;
l = length(x);
for idx = 1:l
	sum = sum + x(idx);
end
end