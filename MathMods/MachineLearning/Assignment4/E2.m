%%	Exercise 2 - Digit classification with LDA


%%	Initialize
clear all; clf; clc;
digits = [2 9];

%	Train
load('usps_train.mat');
trList = find(train_label == digits(1) | train_label == digits(2));
x_train = double(train_data(trList,:));
y_train = double(train_label(trList));

%	Test
load('usps_test.mat');
tsList = find(test_label == digits(1) | test_label == digits(2));
x_test = double(test_data(tsList,:));
y_test = double(test_label(tsList));


%%	Part (a)   	Training

cls = ClassificationDiscriminant.fit(x_train,y_train);
%	Get classes
c1 = cls.Coeffs(1,2).Class1;   c2 = cls.Coeffs(1,2).Class2;
K = cls.Coeffs(1,2).Const;
L = cls.Coeffs(1,2).Linear;

%	Visualize the training example number 155.
dig = reshape(L, 16, 16);
imagesc(dig);
colorbar();  % colormap('gray');


%%	Part (b)   	Predict
f = @(x) (K + x * L);
y_pred_lda = @(x, c1, c2) ((c1 + c2)/2 + (c1 - c2)/2 * sign(f(x)));
lossLDA = loss01(y_pred_lda(x_test, c1, c2), y_test)
loss5NN = loss01(knnClassify(x_train, y_train, x_test, 5), y_test)
