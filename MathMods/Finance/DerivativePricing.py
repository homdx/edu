import math

from os import system    #  TODO: Only to use clear. Remove afterwards.
system('clear')


r = 0.07
s0 = 100.0
T = 1.0

sigma2k090 = {'sigma':0.2, 'k': 90.0}
sigma4k090 = {'sigma':0.4, 'k': 90.0}
sigma2k110 = {'sigma':0.2, 'k':110.0}
sigma4k110 = {'sigma':0.4, 'k':110.0}


def	pos(x):
	return(x if (x > 0) else 0)

def max(x, y):
	return(x if (x > y) else y)


def getTree(s0, sigma, n, t):
	"""Generate the tree of stock prices"""
	
	s = [[0]*(i+1) for i in range(n+1)]
	s[0][0] = s0
	
	dt = t / n
	u = math.exp(sigma * math.sqrt(dt))
	
	for i in range(1, n+1):
		for j in range(i+1):
			s[i][j] = s0 * u ** (-i + 2*j)
	
	return(s)


def callPrices(S0, K, sigma, n, t):
	"""Generate the prices"""
	
	C = [[0]*(i+1) for i in range(n)]
	dt = t / n
	# R = math.exp(r * dt)
	R = 1 + r / 4
	p = (R - d) / (u - d)
	
	# Obtain the binomial tree
	S = getTree(S0, r, u, d, n, t)
	
	# Call prices for the terminal nodes
	for j in range(n):
		C[n-1][j] = pos(S[n-1][j] - K)
	
	for i in range(n-1, -1, -1):
		for j in range(i+1):
			C[i-1][j] = max(S0 * u ** j * d ** (i-j), 1/R * (p * C[i][j+1] + (1-p) * C[i][j]))
	return(C)


def EuropeanAsianCallSP(s0, k, sigma, n, t):
	dt = t / n
	u = math.exp(sigma * math.sqrt(dt))
	d = 1 / u
	R = math.exp(r * dt)
	p = (R - d) / (u - d)
	
	# Obtain the binomial tree
	s = getTree(s0, sigma, n, t)
	
	#  The tree for averages
	sp = [[[]]*(i+1) for i in range(n+1)]
	
	#  Averages for the nodes at maturity (i = n)
	#  ------------------------------------------
	#  Average for N(n,0)
	a = s0 / (n + 1) * (1 - d**(n+1)) / (1 - d)
	sp[n][ 0] = [(a, pos(a - k))]
	#  Average for N(n,n)
	a = s0 / (n + 1) * (1 - u**(n+1)) / (1 - u)
	sp[n][-1] = [(a, pos(a - k))]
	#  Averages for all the other nodes
	for j in range(1, n):
		a_min = s0 / (n+1) * ((1-d**(n-j+1)) / (1-d) + u*d**(n-j) * (1-u**j) / (1-u))
		a_max = s0 / (n+1) * ((1-u**(j+1)) / (1-u) + d*u**j * (1-d**(n-j)) / (1-d))
		if k < a_min:
			sp[n][j] = [(a_min, a_min - k), (a_max, a_max - k)]
		elif k > a_max:
			sp[n][j] = [(a_min, 0), (a_max, 0)]
		else:
			sp[n][j] = [(a_min, 0), (k, 0), (a_max, a_max - k)]
	
	
	#  Averages for the nodes (i,j) forall i != n
	#  ------------------------------------------
	for i in range(n-1, n-2, -1):
		
		#  SP for N(i,0)
#		a = s0 / (i + 1) * (1 - d**(i+1)) / (1 - d)
#		print(i, len(sp[i+1][1]), len(sp[i+1][0]))
#		print(sp[i+1][1][ 0][1], sp[i+1][0][0][1])
#		sp[i][ 0] = [(a, 1/R * (p * sp[i+1][1][ 0][1] + (1-p) * sp[i+1][0][ 0][1]) )]    #  Add max here for American
		
		#  SP for N(i,i)
#		a = s0 / (i + 1) * (1 - u**(i+1)) / (1 - u)
#		sp[i][-1] = [(a, 1/R * (p * sp[i+1][1][-1][1] + (1-p) * sp[i+1][i][-1][1]) )]    #  Add max here for American
		
		
		#  SP for all other nodes N(i,j) s.t. j != 0,i
		for j in range(0, i+1):
			
			#  'B'
			#  Initialize sp_b
			sp_b = list()
			for (a, v_a) in sp[i+1][ j ]:    # TODO: replace by set difference
				b = ((i+2)*a - s0*u**(-i+2*j-1)) / (i+1)
				b_up = ((i+1)*b + s0*u**(-i+2*j+1)) / (i+2)    #  TODO: Error: Check values. Cannot be!!!
				up_points = sp[i+1][j+1]
				if up_points[0][0] < b_up < up_points[-1][0]:
					for kb in range(len(up_points)-1):    #  TODO: replace by binary search
						if up_points[kb][0] < b_up < up_points[kb+1][0]:
							v_b_up = (up_points[kb+1][1] - up_points[kb][1]) / (up_points[kb+1][0] - up_points[kb][0]) * (b_up - up_points[kb][0]) + up_points[kb][1]  #  Value of B_up
							sp_b.append((b_up, 1/R * (p * v_b_up + (1-p) * v_a)))
			
			#  'C'
			#  Initialize sp_c
			sp_c = list()
			for (a, v_a) in sp[i+1][j+1]:
				c = ((i+2)*a - s0*u**(-i+2*j-1)) / (i+1)
				c_dn = ((i+1)*c + s0*d**(-i+2*j+1)) / (i+2)
				dn_points = sp[i+1][j]
				if up_points[0][0] < b_up < up_points[-1][0]:
					for kc in range(len(dn_points)-1):    #  TODO: replace by binary search
						if dn_points[kc][0] < c_dn < dn_points[kc+1][0]:
							v_c_dn = (dn_points[kc+1][1] - dn_points[kc][1]) / (dn_points[kc+1][0] - dn_points[kc][0]) * (c_dn - dn_points[kc][0]) + dn_points[kc][1]  #  Value of C_dn
							sp_c.append((c_dn, 1/R * (p * v_c_dn + (1-p) * v_a)))
			
			#  Aggregate 'B's and 'C's
			sp[i][j] = sorted(sp_b + sp_c, key = lambda sp: sp[0])
			print("sp_b =", len(sp_b), "sp_c =", len(sp_c), "sp[i][j] =", len(sp[i][j]))
			print("sp[i+1][j] =", len(sp[i+1][j]), "sp[i+1][j+1] =", len(sp[i+1][j+1]), "sp[i][j] =", len(sp[i][j]))
	
	return sp


print(getTree(s0 = 100.0, sigma = 0.2, n = 7, t = 1))
eac = EuropeanAsianCallSP(s0 = 100.0, k = 90.0, sigma = 0.2, n = 7, t = 1)
print("Singular Points")
for (i,ptj) in enumerate(eac):
	if i in range(7):
		break
	print("i =", i)
	for j in ptj:
		print(j)
