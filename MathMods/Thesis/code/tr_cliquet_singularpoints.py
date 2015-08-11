__author__ = 'Sudip Sinha'


from math import exp, sqrt, ceil, floor, log


# Market
r = 0.03

# Underlying
# s0 = 50.
sigma = 0.2
q = 0.0

# Derivative
t = 5.
nobs = 5
(f_loc, c_loc, f_glob, c_glob) = (0., 0.08, 0.16, float("inf"))

# k = 60.
# am = False

# Computation
m = 1000
# mach_eps = 65536 * (7 / 3 - 4 / 3 - 1)
# h = 1e-4
# ub = True


# @profile
def sp_cliquet(
		r : float, q : float, sigma : float, t : float,
		m : int, nobs : int,
		f_loc : float, f_glob : float,
		c_loc : float, c_glob : float,
		mach_eps = 65536 * (7/3 - 4/3 - 1)) -> float:
	"""Prices of an Cliquet call option using the singular point method"""

	# Time period for each option
	dt = t / nobs

	# Checks
	f_glob = max( nobs * f_loc, f_glob )
	c_glob = min( nobs * c_loc, c_glob )

	# Singular points for maturity
	sp_i = [(nobs * f_loc, f_glob), (f_glob, f_glob), (c_glob, c_glob), (nobs * c_loc, c_glob)]
	print(sp_i)

	# Unique volatility and rates of interest assumed
	unique_r = True
	unique_v = True
	unique_fc = True

	if unique_r:
		R = exp((r - q) * dt)    # Effective interest rate

	if unique_v:
		# Up and down factors
		u = exp(sigma * sqrt(dt))
		# d = 1. / u

	if unique_r and unique_v:
		# Risk neutral probability, discounted by R
		p_u = ( R * u - 1 ) / ( u * u - 1 )
		p_d = 1. - p_u

		if p_u < 0. or p_d < 0. or p_u > 1. or p_d > 1.:
			raise ValueError

	if unique_fc:
		# Possible range of return
		j_min = floor( log( f_loc + 1 ) / ( 2 * sigma * sqrt(dt)) + m / 2 )
		j_max = ceil ( log( c_loc + 1 ) / ( 2 * sigma * sqrt(dt)) + m / 2 )
		j0 = j_max - j_min
	print('m = {m}, j_min = {j_min}, j_max = {j_max}'.format(m = m, j_min = j_min, j_max = j_max))

	# TODO: Put condition: If all are unique...
	# ret denotes R'
	ret = [0 for i in range( j0 + 1 )]
	ret[0] = f_loc
	ret[j0] = c_loc
	for j in range(1, j0):
		ret[j] = u ** ( -m + 2 * ( j + j_min )) - 1

	# prb denotes p'
	prb = [0 for i in range( j0 + 1 )]
	prb[0] = 0
	for j in range(0, j_min + 1):
		prb[0] += combn( m, j ) * ( p_u ** j ) * ( p_d ** ( m - j ) )
	prb[j0] = 0
	for j in range( j_max, m + 1 ):
		prb[j0] += combn( m, j ) * ( p_u ** j ) * ( p_d ** ( m - j ) )
	for j in range( 1, j0 ):
		prb[j] = combn( m, j ) * ( p_u ** j ) * ( p_d ** ( m - j ) )

	# Go back in time, one step at a time.
	for i in range( nobs - 1, -1, -1 ):
		sp = []

		# Obtain a sorted list for B
		b_list = []
		for l in range( 0, len(sp_i) ):
			for j in range( 0, j0 + 1 ):
				b = sp_i[l][0] - ret[j]
				if ( i * f_loc ) <= b <= ( i * c_loc ):
					close = False
					for rsz in b_list:
						if abs(rsz - b) < mach_eps:
							close = True
							break
					if not close:
						b_list.append(b)
		b_list = sorted(b_list)

		for b in b_list:
			v = 0
			for k in range( 0, j0 + 1 ):
				rszp = b + ret[k]    # Running sum Z corresponding to B at the next time step.
				if rszp <= sp_i[0][0]:
					vk = sp_i[0][1]
				elif rszp >= sp_i[-1][0]:
					vk = sp_i[-1][1]
				else:
					for idx in range( 0, len(sp_i) ):
						if sp_i[idx][0] <= rszp < sp_i[idx + 1][0]:
							vk = sp_i[idx][1] + \
								 ( sp_i[idx + 1][1] - sp_i[idx][1] ) / \
								 ( sp_i[idx + 1][0] - sp_i[idx][0] ) * \
								 ( rszp - sp_i[idx][0] )
				v += prb[k] * vk
			v /= R
			sp.append( (b, v) )

		sp_i = sp

	return(sp_i[0][1])


def factorial(n):
	f = 1
	for i in range( 2, n + 1 ):
		f *= i
	return f


def combn(n: int, r: int):
	if r > ( n / 2 ):
		r = n - r
	c = 1
	for i in range( n, n - r, -1 ):
		c *= i
	return c // factorial(r)


print(sp_cliquet(r=r, q=q, sigma=sigma, t=t, m=m, nobs=nobs, f_loc=f_loc, f_glob=f_glob, c_loc=c_loc, c_glob=c_glob))
