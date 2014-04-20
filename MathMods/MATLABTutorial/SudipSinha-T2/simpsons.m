function [ int ] = simpsons( f, a, b, n )
%simpsons Composite Simpson's rule for integration

h = (b-a)/n;
c = h/2;
x = a:c:b;
y = f(x);
int = (h/6)*(y(1) + 2*sum(y(3:2:2*n-1)) + 4*sum(y(2:2:2*n)) + y(2*n+1));

end
