%%	Exercise 5 - PCA on USPS data


%%	Initialize
clear all;  close all;  clc;


%%	Part (a)
load('usps_train_complete.mat');

% X = double(train_data);
% n = size(X, 1);
% X_centered = X - repmat(mean(X), [n 1]);
% Covariance = X_centered'*X_centered;
% [V D] = eig(Covariance);
coeff = pca(train_data);

pca_v1 = reshape(coeff(:,end),   [16 16]);
pca_v2 = reshape(coeff(:,end-1), [16 16]);

figure(1); 
subplot(1,2,1);  imagesc(pca_v1);  title('Principal Component 1');  colormap(gray);
subplot(1,2,2);  imagesc(pca_v2);  title('Principal Component 2');  colormap(gray);


%%	Part (b)
digit5_labels = find(train_label == 5);
random_array = randperm(length(digit5_labels));
selected_labels = digit5_labels(random_array(1:3));
X = train_data(selected_labels,:);

n = size(X,1);
X_centered = X - repmat(mean(X),n,1);
Covariance = X_centered'*X_centered;
[coeff, ~] = eig(Covariance);

pca_v1 = coeff(:,end);
pca_v2 = coeff(:,end-1);

projection_1 = (pca_v1 * pca_v1') * X_centered' + repmat(mean(X)', [1 3]);
projection_2 = [pca_v1 pca_v2] * [pca_v1 pca_v2]' * X_centered' + repmat(mean(X)', [1 3]);

figure(2);
for i = 1 : 3
    subplot(3, 3, i);
    imagesc(reshape(X(i,:), [16 16]));
    subplot(3, 3, 3 + i);
    imagesc(reshape(projection_1(:,i), [16 16]));
    subplot(3, 3, 6 + i);
    imagesc(reshape(projection_2(:,i), [16 16]));
end
colormap(gray);
