%%	Exercise 4 - Decision Boundary

%%Initialization
clear all;  clc;

n=1000;
p1 = 0.4;
p2 = 0.6;
[X, Y] = mixGaussian2d(n,p1,p2);
n1 = length(find(Y==1));

%%	Part (a)
mean_1 = mean(X(1:n1,:));
mean_2 = mean(X(n1+1:n,:));
variance1 = var(X(1:n1,:));
variance2 = var(X(n1+1:n,:));

%	Plot the points of every class
plot(X(1:n1,1),X(1:n1,2),'b*')
hold on
plot(X(n1+1:n,1),X(n1+1:n,2),'r*')
grid on
xlabel('X_1');  ylabel('X_2');

%%	Part (b)
syms x1;
syms x2;
x = [x1;x2];

% calculate all the vectors and matrices that we need according to Duda
% Book, chapter 2.6.3
s1 = cov(X(1:n1,:));
s2 = cov(X(n1+1:n,:));
m1 = mean_1;
m2 = mean_2;
inverse1 = pinv(s1);
inverse2 = pinv(s2);
W1 = (-1/2)*inverse1;
W2 = (-1/2)*inverse2;
w1 = inverse1*m1';
w2 = inverse2*m2';
w10 = -((1/2)*dot(m1',w1))-(1/2).*log(det(s1))+log(p1);
w20 = -((1/2)*dot(m2',w2))-(1/2).*log(det(s2))+log(p2);

g1 = W1(1,1)*x1^2+(W1(1,2)+W1(2,1))*x1*x2+W1(2,2)*x2^2 + w1'*x + w10;
g2 = W2(1,1)*x1^2+(W2(1,2)+W2(2,1))*x1*x2+W2(2,2)*x2^2 + w2'*x + w20;
g = g1 - g2;
% find the coefficients of the 2 variables polynomial
C = coeffs(g);
V = sym2poly(C); %V has the form [const x2 x2^2 x1 x1*x2 x1^2]
% plot the decision boundary
f = @(x,y) V(6)*x^2 + V(5)*x*y +V(4)*x + V(3)*y^2 + V(2)*y +V(1);
hold on
ezplot(f,[-10, 10, -10, 10])
xlabel('X_1');  ylabel('X_2');
title('Analytic Decision Boundary');
grid on;
