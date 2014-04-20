function [ x, y ] = euler_mod( odefun, x0, y0, xf, n )
%euler_mod Modified Euler method for solving ODEs

h = (xf-x0)/n;
x = [x0 zeros(1,n)];
y = [y0 zeros(1,n)];

for i = 1:n
	x(i+1) = x(i) + h;
	y(i+1) = y(i) + h*odefun(y(i));
	y(i+1) = y(i) + (h/2)*(odefun(y(i)) + odefun(y(i+1)));
end

end
