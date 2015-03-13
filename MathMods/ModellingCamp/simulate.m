function [ meanC ] = simulate( nStickers, nPlayers, packageSize, iterations )

tic;
s = 0;

for i = 1:iterations
	s = s + expectedCost( nStickers, nPlayers, packageSize );
end

meanC = s / iterations;
toc;

end
