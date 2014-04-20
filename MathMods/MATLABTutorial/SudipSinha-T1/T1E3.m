%	Test 1 - Exercise 3
clc;
clear;


%	Define the matrix and vector.
A = [ 1 -3 -5;
	  3 11 -9;
	 -1  1  6];
b = [2 4 5];

%	Backing up 'b' to check correctness later.
bbackup = b;

%	Detecting matrix type.
if A == A'
	fprintf('The matrix A is symmetric.\n')
else
	fprintf('The matrix A is not symmetric.\n')
end
if A == A'
	[n,m]=size(A);
	for i = 1:n
		detAi = det(A(1:i, 1:i));
		if detAi <= 0
			fprintf('The matrix A is not positive definite.\n');
			return
		else
			fprintf('The matrix A is simmetric and positive definite.\n');
		end
	end
end

%	Use the method of Gaussian elimination to solve the problem.
[U, b] = el_gauss(A, b);
x = sol_triang_sup(U, b);
[n1, n2] = size(x);
x'

%	Justification
%	If our solution is correct, we must get A*x' = b'.
if A * x' == bbackup'
	disp('Solution is justified');
end
