function [ root ] = newtonRaphson( fn, x0, tol, itermax )
%NewtonRaphson Find the root of the function 'fn' between 'a' and 'b' using the Newton-Raphson method.
%   'tol' is the tolerance
%	'itermax' is the maximum number of iteration

syms x;

iter = 0;
root = x0;
f = @(x) (fn + x);
f = f(x);
f = @(x) (subs(f, x));
Df = differentiate(f);

while ((iter < itermax) && ((iter == 0) || (abs(root - rootOld) > tol)) && (fn(root) > tol))
%	fprintf('Iteration=%i: root=%f\n', iter, root);
	iter = iter + 1;
	rootOld = root;
	fRoot = fn(root);
%	dfRoot = (fn(root+tol) - fn(root-tol)) / (2 * tol);
	dfRoot = double(Df(root));
	root = root - fRoot / dfRoot;
end

end
