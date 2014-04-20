function [ df ] = differentiate( f )
%differentiate Sybolically differentiates f and returns the derivative.
%   The returned function is an anonymous function.

syms x;
if (~isequal(class(f), 'function_handle'))
	fprintf('ERROR: ''f'' must be a function_handle');
	return;
end

g = diff(f(x));
df = @(x) (subs(g, x));

end
