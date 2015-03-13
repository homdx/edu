import math

def getTree(S0, K, r, sigma, n, t):
	"""Generate the tree of S"""
	
	S = [[0]*(i+1) for i in range(n+1)]
	S[0][0] = S0
	
	dt = t / n
	u = math.exp(sigma * math.sqrt(dt))
	
	for i in range(1, n+1):
		for j in range(i+1):
			S[i][j] = S0 * u ** (-i + 2*j)
	
	return(S)


def callPrices(S0, K, r, u, d, n, t):
	"""Generate the prices"""
	
	C = [[0]*(i+1) for i in range(n)]
	dt = t / n
	# R = math.exp(r * dt)
	R = 1 + r / 4
	p = (R - d) / (u - d)
	
	# Obtain the binomial tree
	S = getTree(S0, K, r, u, d, n, t)
	
	# Call prices for the terminal nodes
	for j in range(n):
		C[n-1][j] = pos(S[n-1][j] - K)
	
	for i in range(n-1, -1, -1):
		for j in range(i+1):
			C[i-1][j] = max(S0 * u ** j * d ** (i-j), 1/R * (p * C[i][j+1] + (1-p) * C[i][j]))
	
	print(S)
	print(C)
	
	return(C)


def	pos(x):
	return(x if (x > 0) else 0)

def max(x, y):
	return(x if (x > y) else y)


print(getTree(S0 = 50, K = 51, r = 0.05, sigma = 0.1, n = 2, t = 1))
