function [ pred ] = knnClassify( train_data, train_label, test_data, k )
%knnClassify The k-nearest neighbour classifier
%	Classifies the points of the test_data according to the
%	k-nearest neighbours algorithm.

%%	Initialization
[mTrain, nTrain] = size(train_data);
[mTest,  nTest] = size(test_data);
pred = zeros(mTest, 1);
dist = zeros(mTrain, 1);
classes = unique(train_label);
nClass = zeros(length(classes), 1);

%%	Body
% Loop over test data
for iTest = 1:mTest

	%	Loop over training data
	for iTrain = 1:mTrain
		dist(iTrain) = norm(train_data(iTrain,:) - test_data(iTest,:));
	end
	
	[~, sIdx] = sort(dist);	% Find out the closest elements.
	sIdx = sIdx(1:k);	% Keep only the min 'k' elements.
	
	% N(Points in Class 'idx').
	for idx = 1:length(classes)
		nClass(idx) = length(find(train_label(sIdx) == classes(idx)));
	end
	
	[~, MIdx] = max(nClass);
	pred(iTest) = classes(MIdx);
end

end
