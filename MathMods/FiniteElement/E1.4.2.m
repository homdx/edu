%%	Group 2 - Assignment 1
%	Anirudh Gupta     -  6607652
%	Maria Panoukidou  -  6607946
%	Sudip Sinha       -  6609779

%%	Initialize
clc; clear all;
syms x;

N = 2^4;
h = 2/N;
xd = 0:h:2;
nn = length(xd);

K = zeros(nn,nn);
u = (64*pi*pi+1)*cos(8*pi*x);
f = zeros(nn);
b = zeros(nn,1);

%%	Calculate K and b
for i = 1:N
	N1 = (1/(xd(i+1)-xd(i)))*(x-xd(i));
	N2 = (1/(xd(i)-xd(i+1)))*(x-xd(i+1));
	N = [N1 N2];
	dN = diff(N,x);
	kef = dN'*dN + N'*N;
	ke = int(kef,0,2);
	temp = zeros(nn, nn);
	temp(i:i+1,i:i+1) = ke;
	K = K + temp;
	f1 = u*N1;
	f2 = u*N2;
	if1 = int(f1,xd(i),xd(i+1));
	if2 = int(f2,xd(i),xd(i+1));
	b(i) = b(i) + double(if1);
	b(i+1) = b(i+1) + double(if2);
end

cond(K)
uh = K \ b;
