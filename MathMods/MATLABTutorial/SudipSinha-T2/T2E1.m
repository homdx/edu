%	Test 2 - Exercise 1
clc;
clear;

n = 100;
x0 = 0;
xf = 2;
y0 = 0.3;
f = @(x, y) (exp(x).*sin(2*x));

[xef, yef] = euler_forward(f, x0, y0, xf, n);
[xem, yem] = euler_mod(f, x0, y0, xf, n);

h = (xf-x0)/n;
xrk(1) = x0;
yrk(1) = y0;

for i = 1:n
	k1(i) = f(xrk(i),       yrk(i));
	k2(i) = f(xrk(i) + h/2, yrk(i) + h/2 * k1(i));
	k3(i) = f(xrk(i) + h/2, yrk(i) + h/2 * k2(i));
	k4(i) = f(xrk(i) + h,   yrk(i) + h   * k3(i));
	
	xrk(i+1) = xrk(i) + h;
	yrk(i+1) = yrk(i) + h/6 * (k1(i) + 2 * k2(i) + 2 * k3(i) + k4(i));	
end

plot(xef,yef,'b', xem,yem,'g', xrk,yrk,'r');
legend('Euler (Forward)', 'Euler (Modified)', 'RK2'); grid on;
xlabel('x'); ylabel('y'); title('Solution of Cauchy problem');
