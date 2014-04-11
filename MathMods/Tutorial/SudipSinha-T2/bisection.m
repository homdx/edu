function [ root ] = bisection( fn, a, b, tol, itermax )
%bisection Find the root of the function 'fn' between 'a' and 'b' using the bisection method.
%   'tol' is the tolerance
%	'itermax' is the maximum number of iteration

if (nargin < 3)
	disp('Error: you must mention the function and the range of data.');
	return;
end

%	Preliminary checks.
if ((fn(a) - 0) < tol)
	root = a;
elseif ((fn(b) - 0) < tol)
	root = b;
end

iter = 0;
while ((abs(b-a) > tol) && (iter < itermax))
	iter = iter + 1;
	c = (a + b) / 2;
	if (abs(fn(c)) < tol || ((b - a) < tol))
		root = c;
		break;
	elseif (fn(a) * fn(c) < 0)
		b = c;
	else
		a = c;
	end
	root = c;
end
end
