function [ pred ] = knnClassify( train_data, train_label, test_data, k )
%knnClassify The k-nearest neighbour classifier
%	Classifies the points of the test_data according to the
%	k-nearest neighbours algorithm.

mTrain = size(train_data, 1);
mTest  = size(test_data,  1);

% size(dist) = mTrain x mTest
xx = sum(train_data .^ 2, 2);
yy = sum(test_data  .^ 2, 2);
xy = train_data * test_data';
dist =  repmat(xx, [1  mTest]) + repmat(yy', [mTrain 1]) - 2*xy;

[~, sIdx] = sort(dist, 1);	% Find out the closest elements

train_label = repmat(train_label, [1 mTest]);

pred = mode(train_label(sIdx(1:k, :)), 1);
pred = pred';

end
