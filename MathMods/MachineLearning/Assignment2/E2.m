%%	Exercise 2 - Bayes classifier

%%	Initialize
clear all; clc;
[X, Y] = mixGaussian1d(1000, 0.5, 0.5, 0, 6, 1, 2);
% [X, Y] = mixGaussian1d(1000, 0.1, 0.6, 0, 6, 1, 2);

%% Part (a)

%	Description of the data.
%	According to p1 and p2, it generates normally distributed data bases on
%	mu and sigma values for two classes.

figure(1);  clf;

%	Scatter
subplot(2, 2, 1);
scatter(X, Y, '.');
grid on;  title('Scatter');
xlabel('X');  ylabel('Y');  ylim([0.9 2.1]);

%	Marginal
subplot(2, 2, 2);
hold all;  grid on;  title('Marginal');
[countC, binsX] = hist(X, 30);
PX = countC / size(X, 1);
plot(binsX, PX, 'r.-');
xlabel('X');  ylabel('P(X)');

%	Class conditionals.
subplot(2, 2, 3);
hold all;  grid on;  title('Class Conditionals');
%	P(X | Y = 1)
XY1 = X(Y == 1);
[countC, binsX] = hist(XY1, binsX);
PXY1 = countC / size(XY1, 1);
plot(binsX, PXY1, 'ro-');
%	P(X | Y = 2)
XY2 = X(Y == 2);
[countC, binsX] = hist(XY2, binsX);
PXY2 = countC / size(XY2, 1);
plot(binsX, PXY2, 'b*-');
xlabel('X');  ylabel('P(X)');  legend('P(X | Y = 1)', 'P(X | Y = 2)');

%	Priors
PY1 = length(find(Y == 1)) / length(Y);
PY2 = length(find(Y == 2)) / length(Y);
fprintf('P(Y = 1) = %f, P(Y = 2) = %f\n', PY1, PY2);

%	Posteriors
PY1X = (PXY1 .* PY1) ./ PX;
PY2X = (PXY2 .* PY2) ./ PX;

subplot(2, 2, 4);
hold all;  grid on;  title('Posteriors');
%	P(Y = 1 | X)
plot(binsX, PY1X, 'ro-');
%	P(Y = 2 | X)
plot(binsX, PY2X, 'b*-');
xlabel('X');  ylabel('P(X)');  legend('P(Y = 1 | X)', 'P(Y = 2 | X)');



%%	Part (b)
Xx = [0 1 2 2.5 3 4.5];
PXx = interp1(binsX, PX, Xx);	% From the marginals graph.

%	MLE
PXxY1 = interp1(binsX, PXY1, Xx);	% From the class conditionals graph.
PXxY2 = interp1(binsX, PXY2, Xx);	% From the class conditionals graph.
mle = ones(1, length(Xx));
mle(PXxY2 > PXxY1) = 2;
mle

%	Bayes
PY1Xx = interp1(binsX, PY1X, Xx);	% From the posteriors graph.
PY2Xx = interp1(binsX, PY2X, Xx);	% From the posteriors graph.

bay = ones(1, length(Xx));
bay(PY2Xx > PY1Xx) = 2;
bay
