%	Test 1 - Exercise 2
clf;
clc;
clear;


%	Figure 1
figure(1);

subplot(1, 3, 1);
x = -2:0.01:2; y = x.^2 + 2*x + 1;
plot(x, y);
title('Graph of y = x^2 + 2x + 1');
xlabel('x'); ylabel('y');
grid on;

subplot(1, 3, 2);
x = 0:0.01:5; y = sin(x) .* exp(x);
plot(x, y);
title('Graph of y = e^xsin(x)');
xlabel('x'); ylabel('y');
grid on;

subplot(1, 3, 3);
x = -10:0.01:10; y = cos(x) .* x;
plot(x, y);
title('Graph of y = x.cos(x)');
xlabel('x'); ylabel('y');
grid on;


%	Figure 2
figure(2);

subplot(1, 3, 1);
[x, y] = meshgrid(-2:0.01:2); z = x.^2+y.^2;
mesh(x, y, z);
title('z = x^2 + y^2');
xlabel('x'); ylabel('y'); zlabel('z');
grid on; colorbar;

subplot(1, 3, 2);
x = -5:0.01:5; y = -6:0.01:6;
[x, y] = meshgrid(x, y); z = x.^2 - y.^2;
mesh(x, y, z);
title('z = x^2 - y^2');
xlabel('x'); ylabel('y'); zlabel('z');
grid on; colorbar;

subplot(1, 3, 3);
[x, y] = meshgrid(-10:0.3:10); z = sin(sqrt(x.^2+y.^2))/sqrt(x.^2+y.^2);
mesh(x, y, z);
title('z = sin(sqrt(x^2+y^2))/sqrt(x^2+y^2)');
xlabel('x'); ylabel('y'); zlabel('z');
grid on; colorbar;


%	Isolines
figure(3);

subplot(1, 3, 1);
[x, y] = meshgrid(-2:0.01:2); z = x.^2+y.^2;
contour(x, y, z);
title('z = x^2 + y^2');
xlabel('x'); ylabel('y');
grid on; colorbar;

subplot(1, 3, 2);
x = -5:0.01:5; y = -6:0.01:6;
[x, y] = meshgrid(x, y); z = x.^2 - y.^2;
contour(x, y, z);
title('z = x^2 - y^2');
xlabel('x'); ylabel('y');
grid on; colorbar;

subplot(1, 3, 3);
[x, y] = meshgrid(-10:0.3:10); z = sin(sqrt(x.^2+y.^2))/sqrt(x.^2+y.^2);
contour(x, y, z);
title('z = sin(sqrt(x^2+y^2))/sqrt(x^2+y^2)');
xlabel('x'); ylabel('y');
grid on; colorbar;
