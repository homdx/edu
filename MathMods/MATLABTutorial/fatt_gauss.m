function [L,U] = fatt_gauss(A)
	[nrow, ncol] = size(A);
	L = eye(nrow, ncol);
	for k = 1:(nrow-1)
		for i = (k+1):nrow
			L(i,k) = A(i,k) / A(k,k);
			for j = k:ncol
				A(i,j) = A(i,j) - L(i,k) * A(k,j);
			end
		end
	end
	U = A;
end