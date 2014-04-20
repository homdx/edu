%%	Group 2 - Assignment 1
%	Anirudh Gupta     -  6607652
%	Maria Panoukidou  -  6607946
%	Sudip Sinha       -  6609779

%%	Initialize
clc; clear all;
syms t;

N = 16;
ll = 0; ul = 2;	% Lower limit and upper limit of Omega.
A = zeros(N+1);
b = zeros(N+1,1);

%% Calculate A and b
for i = 0:N
	phi_i = zeros(1, N+1);
	phi_i(N+1-i) = 1;
	for j = 0:N
		phi_j = zeros(1, N+1);
		phi_j(N+1-j) = 1;
		a1 = conv(polyder(phi_i), polyder(phi_j));
		a2 = conv(phi_i, phi_j);
		diff = length(a2) - length(a1);
		if (diff > 0)
			a1(1+diff:end+diff) = a1;
			a1(1:diff) = zeros(1, diff);
		end
		a = a1 + a2;
		I = polyint(a);
		A(i+1,j+1) = polyval(I,ul) - polyval(I,ll);
	end
end

fprintf('Condn(A) = %f\n', cond(A));


syms x;
f = 64*pi*pi*cos(8*pi*x)+cos(8*pi*x);
b = zeros(N+1, 1);
for i = 0:N
	b(i+1) = int(f * x^i, x, 0, 2);
end
b = double(b);

uhi = A \ b;

%%	Calculate the interpolating function
uh = 0;
for i = 0:N
	uh = uh + uhi(i+1) * x^i;
end

xi = 0:0.1:2;

%%	Plot
hold all;
p1 = ezplot(uh, [0,2]);
set(p1,'Color','red');
p2 = ezplot(cos(8*pi*x), [0,2]);
set(p2,'Color','blue');
legend('u_h', 'cos(8\pix)');
grid on;
