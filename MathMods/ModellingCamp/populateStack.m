function [ stack ] = populateStack( nStickers, players, packageSize, stack )
%populateStack    Populate the stack by randomly drawing stickers
%   Detailed explanation goes here

%	Populate the stack
for i = players
	%	Draw 'packageSize' number of stickers randomly.
	draw = 1 + floor(nStickers * rand(packageSize, 1));
	
	%	Fill the stack
	stack(i, draw) = stack(i, draw) + 1;
end

end
