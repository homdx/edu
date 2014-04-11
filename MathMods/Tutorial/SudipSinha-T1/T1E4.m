%	Test 1 - Exercise 4
clf;
clc;
clear;


x = -2:0.01:16;
y = 10 * (1 - exp(-x/4));

figure(1);
plot(x, y, x, 9.8);
title('Graph of 10(1-e^{-x/4})');
xlabel('x'); ylabel('y');
legend('y = 10(1-e^{-x/4})', 'y = 9.8');
grid on;

idx = length(x);
while y(idx) > 9.8+eps
	idx = idx-1;
end

fprintf('Value of x_i when y(x_i)=9.8 is %f\n', x(idx));
%	This can also be verified from the attached figure (T1E4.fig).
