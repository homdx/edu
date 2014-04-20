%	Test 2 - Exercise 3
clc;
clear;

tol = 1e-4;
f = @(x) (exp(x) .* sin(x) + x.^4);
a = 0; b = pi;

hsi = 2 * (((180*tol)/(10*0.375))^(0.25));
htr = sqrt((12*tol)/(10*2));

ntr = round((b-a)/htr);
nsi = round((b-a)/hsi);

itr = trapezoidal(f, a, b, ntr);
isi = simpsons(f, a, b, nsi);

% Analytic calculation of integral.
iAnalyic = 0.1 * (5 + 5*exp(pi) + 2*pi^5);

etr = abs(iAnalyic - itr);
esi = abs(iAnalyic - isi);

fprintf('Function: f(x) = exp(x)sin(x)+x^4\n');
fprintf('Actual value of the integral:\t%f\n', iAnalyic);
fprintf('Method\t\t\tIntegral\t\tError\n');
fprintf('Trapeziod\t\t%f\t\t%f\n', itr, etr);
fprintf('Simpsons\t\t%f\t\t%f\n', isi, esi);

if (etr < esi)
	disp('Trapezoid method is more accurate in this case.');
else
	disp('Simpson''s method is more accurate in this case.');
end
