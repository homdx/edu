n = 5000;	% number of simulations for product level.
level = 25:50;	% product level
cost = 30000 + 2000 * level;

for k = level
	cum_profit = 0;	% total profit for level production
	for m = 1:n
		demand = floor(rand * (50-25) + 25);
		if demand >= k
			income = 4000 * k;
		else
			income = 4000 * demand + 1000 * (k - demand);
		end
		profit = income - cost(k-24);
		cum_profit = cum_profit + profit;
	end
	expected_profit = cum_profit / n;
	p(k-24,1) = k;
	p(k-24,2) = expected_profit;
end

plot(p(:,1), p(:,2), 'o', p(:,1), p(:,2), '-'); grid on;
xlabel('Unit'); ylabel('Profit');
