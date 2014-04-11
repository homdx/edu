%	Test 2 - Exercise 2
clc;
clear;

f = @(x) (x.^4 + 10*x.^3 - 35*x.^2 - 50*x -100);
Df = differentiate(f);
D2f = differentiate(Df);

x = -10:0.1:10;
xi =   1:0.1:10;
%	Part 1
plot(x, f(x), x, Df(x), x, D2f(x)); grid on;
xlabel('x'); ylabel('y'); title('y = f(x)');
legend('f(x)', 'f''(x)', 'f''''(x)');

%	Part 2
cond = true;
%	Hypothesis 1: f(a).f(b) < 0
if (f(1)*f(10) >= 0)
	cond = false;
	disp('Hypothesis 1 is not satisfied.');
end
%	Hypothesis 2: f'(x) ~= 0
if (length(find(Df(xi) == 0)) ~= 0)
	cond = false;
	disp('Hypothesis 2 is not satisfied.');
end
%	Hypothesis 3: f''(x) ? 0 ? f''(x) ? 0
if (length(find(D2f(xi) >= 0)) ~= 0) && (length(find(D2f(xi) >= 0)) ~= length(xi))
	cond = false;
	disp('Hypothesis 3 is not satisfied.');
end
if (cond)
	disp('All hypotheses for the Newton-Raphson Method are satisfied for 1:10.');
else
	disp('Some hypotheses for the Newton-Raphson Method are not satisfied.');
end

%	Part 3
%	From the plot, we see that there is a root in near x=4.
%	Using Newton's method with x0 = 4.
format long;
rNew = newtonRaphson(f, 4, 1e-15, 100)
rbis = bisection(f, -10, 10, 1e-15, 100)

%	COMPARISON: Both the roots come out to be equal.
%	But for the Newton's method, we	need a starting point.
format short;