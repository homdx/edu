function [ int ] = trapezoidal( f, a, b, n )
%trapezoidal Composite trapezoidal rule for integration

h = (b-a)/n;
x = a:h:b;
y = f(x);
int = h*(0.5*y(1) + sum(y(2:n)) + 0.5*y(n+1));

end
