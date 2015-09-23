__author__ = 'Sudip Sinha'


from math import exp, sqrt, ceil, floor, log, factorial


# Market
r = 0.03

# Underlying
# s0 = 50.
sigma = [ ( 0.05 + 0.04 * i ) for i in range( 1, 9 ) ]
q = 0.0

# Derivative
t = 2.
N = 8
(f_loc, c_loc, f_glob, c_glob) = (0., 0.08, 0.16, float("inf"))

# k = 60.
# am = False

# Computation
m = 5
# h = 1e-4


# @profile
def sp_cliquet(
		r : float, q : float, sigma : float, t : float,
		m : int, N : int,
		f_loc : float, f_glob : float,
		c_loc : float, c_glob : float,
		h : float = 0., mach_eps = 65536 * ( 7/3 - 4/3 - 1 ) ) -> float:
	"""Prices of an Cliquet call option using the singular point method"""

	err_p = "ERROR: p_u must be a probability! " + \
	          "We must have dt < ( sigma ** 2 / ( r - q ) ** 2 ) for this. " + \
	          "You may want to decrease the time interval."

	# Time period for each option
	n = m * N
	dt = t / n

	# Checks
	f_glob = max( N * f_loc, f_glob )
	c_glob = min( N * c_loc, c_glob )

	# Singular points for maturity (go to the future)
	sp_i = [(N * f_loc, f_glob),
	        (f_glob, f_glob),
	        (c_glob, c_glob),
	        (N * c_loc, c_glob)]
	if f_glob == N * f_loc:
		sp_i = sp_i[1:]
	if c_glob == N * c_loc:
		sp_i = sp_i[:-1]
	# print(sp_i)

	# Uniqueness of sigma
	if type(sigma) == float:
		unique_v = True
		assert dt < ( sigma ** 2 / ( r - q ) ** 2 ), err_p
	elif type(sigma) == list and len(sigma) == N:
		unique_v = False
		sigma_list = sigma
	else:
		raise ValueError

	# Uniqueness of r
	if type(r) == float:
		unique_r = True
	elif type(r) == list and len(r) == N:
		unique_r = False
		r_list = r
	else:
		raise ValueError

	# Unique f_loc and c_loc
	unique_fc = True

	if unique_r:
		# Effective interest rate for observable time periods
		R_N = exp( ( r - q ) * t / N )
		# Effective interest rate for computational time periods
		R_n = exp( ( r - q ) * t / n )

	if unique_v:
		u = exp( sigma * sqrt( dt ) )    # Up factor
		# d = 1. / u

	if unique_r and unique_v:
		# Risk neutral probability, discounted by R
		p_u = ( R_n * u - 1 ) / ( u * u - 1 )
		p_d = 1. - p_u
		# print('p_u = {p_u}'.format(p_u = p_u))
		assert 0. <= p_u <= 1. and 0. <=  p_d <= 1., \
			"p_u and p_d must be valid probabilities."

	if unique_v and unique_r and unique_fc:
		# Possible range of return
		j_min = max( 0, floor( log( f_loc + 1 ) / ( 2 * sigma * sqrt(dt)) + m / 2 ) )
		j_max = min( m, ceil ( log( c_loc + 1 ) / ( 2 * sigma * sqrt(dt)) + m / 2 ) )
		j0 = j_max - j_min

		# prb denotes p'
		prb = [ 0 for j in range( j0 + 1 ) ]
		prb[0] = 0
		for j in range( 0, j_min + 1 ):
			prb[0] += choose( m, j ) * ( p_u ** j ) * ( p_d ** ( m - j ) )
		prb[j0] = 0
		for j in range( j_max, m + 1 ):
			prb[j0] += choose( m, j ) * ( p_u ** j ) * ( p_d ** ( m - j ) )
		for j in range( j_min + 1, j_max ):
			prb[j - j_min] = choose( m, j ) * ( p_u ** j ) * ( p_d ** ( m - j ) )
		# print('sum(prb) = {p}'.format(p = sum(prb)))
		assert abs(sum(prb) - 1.) < mach_eps, "Sum of PMF must be unity."

		# ret denotes R'
		ret = [ 0 for j in range( j0 + 1 ) ]
		ret[0] = f_loc
		ret[j0] = c_loc
		for j in range( j_min + 1, j_max ):
			ret[j - j_min] = u ** ( -m + 2 * j ) - 1

	# Come back from the future, one step at a time.
	for i in range( N - 1, -1, -1 ):
		sp = []

		# Checks
		if not unique_v or unique_r:
			if not unique_r:
				r = r_list[i]
				# Effective interest rate for observable time periods
				R_N = exp( ( r - q ) * t / N )
				# Effective interest rate for computational time periods
				R_n = exp( ( r - q ) * t / n )
			if not unique_v:
				sigma = sigma_list[i]
				# Up factor
				u = exp( sigma * sqrt( dt ) )
				# Possible range of return
				j_min = max( 0, floor( log( f_loc + 1 ) / ( 2 * sigma * sqrt(dt)) + m / 2 ) )
				j_max = min( m, ceil ( log( c_loc + 1 ) / ( 2 * sigma * sqrt(dt)) + m / 2 ) )
				j0 = j_max - j_min
			assert dt < ( sigma ** 2 / ( r - q ) ** 2 ), err_p

			# Risk neutral probability, discounted by R
			p_u = ( R_n * u - 1 ) / ( u * u - 1 )
			p_d = 1. - p_u
			# print('p_u = {p_u}'.format(p_u = p_u))
			assert 0. <= p_u <= 1., "p_u and p_d must be valid probabilities."

			# prb denotes p'
			prb = [ 0 for j in range( j0 + 1 ) ]
			prb[0] = 0
			for j in range( 0, j_min + 1 ):
				prb[0] += choose( m, j ) * ( p_u ** j ) * ( p_d ** ( m - j ) )
			prb[j0] = 0
			for j in range( j_max, m + 1 ):
				prb[j0] += choose( m, j ) * ( p_u ** j ) * ( p_d ** ( m - j ) )
			for j in range( j_min + 1, j_max ):
				prb[j - j_min] = choose( m, j ) * ( p_u ** j ) * ( p_d ** ( m - j ) )
			# print('sum(prb) = {p}'.format(p = sum(prb)))
			assert abs(sum(prb) - 1.) < mach_eps, "Sum of PMF must be unity."

			# ret denotes R'
			ret = [ 0 for j in range( j0 + 1 ) ]
			ret[0] = f_loc
			ret[j0] = c_loc
			for j in range( j_min + 1, j_max ):
				ret[j - j_min] = u ** ( -m + 2 * j ) - 1

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
				# Running sum Z corresponding to B at the next time step.
				rszp = b + ret[k]
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
			v /= R_N
			sp.append( (b, v) )

		# Approximations
		if i % n == 0 and h > 0.:
			l = 0    # l is the starting index
			while l < len(sp):
				approx = True
				for j in range( len(sp), l + 1, -1 ):
					m = ( sp[j][1] - sp[l][1] ) / ( sp[j][0] - sp[l][0] )
					for k in range( l, j ):
						approx = True
						delta = m * ( sp[k][0] - sp[l][0] ) + sp[l][0] - sp[k][1]
						if delta >= h:
							approx = False
							break
					if approx:
						break
				if approx:
					del sp[l+1:j]
				l += 1

		sp_i = sp

	return sp_i[0][1]


# def choose(n: int, r: int):
# 	assert ( type(n) == int and type(r) == int and 0 <= r <= n ), "'n' and 'r' must be non-negative integers"
#
# 	if r > ( n / 2 ):
# 		r = n - r
# 	c = 1
# 	for i in range( r ):
# 		c *= ( n - i ) / ( r - i )
# 	return c


def choose(n: int, r: int):
	assert ( type(n) == int and type(r) == int and 0 <= r <= n ), \
		"'n' and 'r' must be non-negative integers"

	if r > ( n / 2 ):
		r = n - r
	c = 1
	for i in range( n, n - r, -1 ):
		c *= i
	return c // factorial(r)


print(sp_cliquet(r=r, q=q, sigma=sigma, t=t, m=m, N=N,
                 f_loc=f_loc, f_glob=f_glob, c_loc=c_loc, c_glob=c_glob, h = 1e-6))
