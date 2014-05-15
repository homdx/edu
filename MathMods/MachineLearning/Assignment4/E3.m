%%	Exercise 3 - Multiclass classification with LDA


%%	Initialize
clear all; clc;

classes = 1:4;

load('usps_train_complete.mat');  load('usps_test_complete.mat');
trList = zeros(0, 1);  tsList = zeros(0, 1);
for class = classes
	trList = [trList; find(train_label == class)];
	tsList = [tsList; find(test_label  == class)];
end
%	Train
x_train = double(train_data(trList,:));
y_train = double(train_label(trList));
%	Test
x_test = double(test_data(tsList,:));
y_test = double(test_label(tsList));


%%	Part (a) - One vs All

score = zeros(1 + size(x_train, 2), length(classes));

for i = 1:length(classes)
	y = y_train;  % Temporary variable
	y(y ~= classes(i)) = intmax;
	lda = ClassificationDiscriminant.fit(x_train, y);
	c1 = lda.Coeffs(1,2).Class1;
	if (c1 == classes(i))
		K = lda.Coeffs(1,2).Const;
		L = lda.Coeffs(1,2).Linear;
	else
		K = lda.Coeffs(2,1).Const;
		L = lda.Coeffs(2,1).Linear;
	end
	score(:,i) = [K; L];
end

[~, y_test_pred] = max(add1(x_test) * score, [], 2);
y_test_pred = classes(y_test_pred);
y_test_pred = y_test_pred';
loss01(y_test_pred, y_test)


%%	Part (b) - One vs One

f = @(x, K, L) (K + x * L);

score = zeros(size(x_test, 1), length(classes));

for i = 1:length(classes)
	for j = (i+1):length(classes)
		idx = find(y_train == i | y_train == j);
		lda = ClassificationDiscriminant.fit(x_train(idx,:), y_train(idx));
		c1 = lda.Coeffs(1,2).Class1;  c2 = lda.Coeffs(1,2).Class2;
		K = lda.Coeffs(1,2).Const;  L = lda.Coeffs(1,2).Linear;
		if (c1 == i)
			score(:, i) = score(:, i) + (f(x_test, K, L) > 0);
			score(:, j) = score(:, j) + (f(x_test, K, L) < 0);
		else
			score(:, j) = score(:, i) + (f(x_test, K, L) > 0);
			score(:, i) = score(:, j) + (f(x_test, K, L) < 0);
		end
	end
end

[~, y_test_pred] = max(score, [], 2);
y_test_pred = classes(y_test_pred);
y_test_pred = y_test_pred';
loss01(y_test_pred, y_test)
