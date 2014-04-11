function zn = nextZ(z, n)
	zn = 2^((n-1)/2) * sqrt(1 - sqrt(1 - 4^(1-n)*z^2));
end