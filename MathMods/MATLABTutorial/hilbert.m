n = 10;
KA = zeros(1,10);
err = zeros(1,10);
nn = 1:10;
for n = 1:10
	A = mat_hilbert(n);
	KA(n) = cond(A);
	xex = ones(n,1);
	b = A * xex;
	x = A \ b;
	err(n) = norm(x - xex) / norm(xex);
end
semilogy(nn, KA, 'r', nn(2:10), err(2:10), 'b');
grid on;
legend('KA','err');