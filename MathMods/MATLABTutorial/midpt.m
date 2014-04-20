function [ int ] = midpt( f, a, b, n )
%midpt Composite midpoint rule for integration

h = (b-a)/n;
x = (a+h/2):h:b;
y = f(x);
int = h*(sum(y));

end
