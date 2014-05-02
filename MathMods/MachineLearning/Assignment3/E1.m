%%	Exercise 1 - Linear mapping


%%	Initialization
clear all; clc;
load('Adot.mat');

theta = pi/3;
V = [cos(theta) -sin(theta);sin(theta) cos(theta)];


%%	Part (a)
Y = V * X;

%	Plot X and Y
figure(1);  clf;  hold all;
title('Linear maps V and V^T');

scatter(X(1,:), X(2,:), 'b.');
scatter(Y(1,:), Y(2,:), 'r.');
legend('X', 'Y');

%	The linear mapping V rotates the operand by an angle theta.


%%	Part (b)
Z = V' * Y;
scatter(Z(1,:), Z(2,:), 'm.');
legend('X', 'Y', 'Z');
%	The linear transform V'V is an identity transform.
%	Thus it replaces the 'X' by 'Z'.


%%	Part (c)
D1 = [2 0; 0 2];
D2 = [2 0; 0 1];

W1 = D1 * X;
W2 = D2 * X;

%	Plot X and its transformations
figure(2);  clf;  hold all;
title('Linear maps D_1 and D_2');

scatter(X(1,:) , X(2,:) , 'b.');
scatter(W1(1,:), W1(2,:), 'r.');
scatter(W2(1,:), W2(2,:), 'm.');

legend('X', 'D_1X', 'D_2X');

%	D1 doubles the distance of the point from the origin.
%	D2 doubles only the first co-ordinate and retains the second.


%%	Part (d)
U2 = V' * D2 * V * X;

%	Plot X and its transformations
figure(3);  clf;  hold all;
title('Linear map V^TD_2VX');

scatter(X(1,:) , X(2,:) , 'b.');
scatter(U2(1,:), U2(2,:), 'r.');

legend('X', 'V^TD_2V');

%	V'*D2*V rotates the data points by theta, scales the first coordinate
%	by 2, and rotates back the scaled points.
