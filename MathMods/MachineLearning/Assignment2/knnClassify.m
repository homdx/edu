function [ pred ] = knnClassify( train_data, train_label, test_data, k )
%knnClassify The k-nearest neighbour classifier
%	Classifies the points of the test_data according to the
%	k-nearest neighbours algorithm.

%%	Initialization
mTrain = size(train_data, 1);
mTest  = size(test_data,  1);
pred = zeros(mTest, 1, 'uint8');	% Assuming there are a maximum of 256 classes.
dist = zeros(mTrain, 1);
classes = unique(train_label);
nClass = zeros(length(classes), 1);

%%	Body
% Loop over test data
for iTest = 1:mTest

	%	Loop over training data
	dist = sum((train_data - repmat(test_data(iTest,:), [mTrain 1])) .^ 2, 2);
	
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
