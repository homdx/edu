%	Plotting the sinusoidal functions in one graph.
x = -2*pi : pi/10 : 2*pi;
plot(x, sin(x), 'r*-', x, cos(x), 'b+--')
grid on
legend('sin', 'cos')
xlabel('\theta (rad)')
ylabel('f^n')
title('Sinusoidal functions')

%	Plotting z = y + e^x
[x y] = meshgrid(-2:0.1:2);
z = y + exp(x);
surf(x, y, z);
grid on;
colorbar;
xlabel('x'); ylabel('y'); zlabel('z');
title('Graph of z = y + e^x')
