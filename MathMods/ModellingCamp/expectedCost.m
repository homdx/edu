function [ cost ] = expectedCost( nStickers, nPlayers, packageSize )
%ExpectedCost Summary of this function goes here
%   Detailed explanation goes here

%	Intitialize the stack
stack = zeros([nPlayers, nStickers]);

%	This vector gives the users who have still not completed the game
playing = 1:nPlayers;

while (sum(playing) > 0)
	%	Draw new stickers
	stack = populateStack( nStickers, playing, packageSize, stack );
	
	%	Swap
	
	%	Check who all have finished
	playing = find(sum(stack ~= 0, 2) ~= nStickers);
end

cost = sum(stack(:)) / nPlayers;

end
