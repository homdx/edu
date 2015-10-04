# ToDo: Low memory implementation
# ToDo: Use assert

__author__ = 'Sudip Sinha'

from math import exp, sqrt, isclose


# @profile
def sp_asian_call( r: float,  # Market
                   s0: float, sigma: float, q: float,  # Underlying
                   k: float, t: float, am: bool = True,  # Derivative
                   n: int = 25, h: float = 0., ub: bool = True  # Computation
                   ) -> list:
	"""Prices of an Asian call option using the singular point method"""

	# Time period
	dt = t / n
	# Effective interest rate
	R = exp( (r - q) * dt )
	
	# Up and down factors
	u = exp( sigma * sqrt( dt ) )
	d = 1. / u
	
	# Risk neutral probability, discounted by R
	p_u = (R * u - 1) / (u * u - 1) / R
	p_d = 1. / R - p_u
	assert p_u >= 0. and R * p_u <= 1., "Probability!"

	# Singular points for nodes N(i+1,*)

	# Start from maturity
	now = [[]] * (n + 1)

	# Singular points for N(n,0)
	a = s0 / (n + 1) * (1 - d ** (n + 1)) / (1 - d)
	now[0] = [(a, max( a - k, 0. ))]

	# Singular points for N(n,n)
	a = s0 / (n + 1) * (1 - u ** (n + 1)) / (1 - u)
	now[-1] = [(a, max( a - k, 0. ))]

	# Singular points for N(n,j)
	for j in range( 1, n ):
		a_min = s0 / (n + 1) * ((1 - d ** (n - j + 1)) / (1 - d)
		                        + d ** (n - j - 1) * (1 - u ** j) /
		                        (1 - u))
		a_max = s0 / (n + 1) * ((1 - u ** (j + 1)) / (1 - u)
		                        + u ** (j - 1) * (1 - d ** (n - j)) /
		                        (1 - d))
		if k < a_min:
			now[j] = [(a_min, a_min - k), (a_max, a_max - k)]
		elif k > a_max:
			now[j] = [(a_min, 0), (a_max, 0)]
		else:
			now[j] = [(a_min, 0), (k, 0), (a_max, a_max - k)]
	display = '\n' + '-' * 64 + '\n'
	for (j, tj) in enumerate(now):
		display = 'N({i},{j}): {tj}\n'.format(i=n, j=j, tj=tj) + display
	display = '-' * 64 + '\n' + display

	
	# Singular points for N(i,j) i != n
	for i in range( n - 1, -1, -1 ):
		future = now
		now = [[]] * (i + 1)

		for j in range( i + 1 ):

			a_min = s0 / (i + 1) * ((1 - d ** (i - j + 1)) / (1 - d)
			                        + d ** (i - j - 1) * (1 - u ** j) / (1 - u))

			a_max = s0 / (i + 1) * ((1 - u ** (j + 1)) / (1 - u)
			                        + u ** (j - 1) * (1 - d ** (i - j)) / (
			                        1 - d))
			# print('a_min = {}, a_max = {}'.format(a_min, a_max))

			if a_min > a_max:
				(a_min, a_max) = (a_max, a_min)  # Jugaad

			# 'B'
			sp_b = list( )  # Initialize sp_b
			for (a, v_a) in future[j]:
				b = ((i + 2) * a - s0 * u ** (-i + 2 * j - 1)) / (i + 1)

				assert len( future[j + 1] ) > 0, "len( future[j+1] ) == 0"

				# Get b_up for b in [a_min, a_max]
				if a_min <= b <= a_max or isclose( b, a_min ) or isclose( b, a_max ):
					b_up = ((i + 1) * b + s0 * u ** (-i + 2 * j + 1)) / (i + 2)
					future_node = future[j + 1]

					v_b_up = 0
					if isclose( b_up, future_node[-1][0] ):
						v_b_up = future_node[-1][1]
					elif isclose( b_up, future_node[0][0] ):
						v_b_up = future_node[0][1]
					else:
						for kb in range( len( future_node ) - 1 ):
							if future_node[kb][0] <= b_up < future_node[kb + 1][0]:
								v_b_up = ((future_node[kb + 1][1] - future_node[kb][1]) /
								          (future_node[kb + 1][0] - future_node[kb][0]) *
								          (b_up - future_node[kb][0]) + future_node[kb][1])
								break
					sp_b.append( (b, p_u * v_b_up + p_d * v_a) )

			# 'C'
			sp_c = list( )  # Initialize sp_c
			for (a, v_a) in future[j + 1]:
				c = ((i + 2) * a - s0 * u ** (-i + 2 * j + 1)) / (i + 1)

				assert len( future[j] ) > 0, "len( future[j] ) == 0"

				# # Uniqueness: Verify that 'C' is not in the set of 'B'
				# if any((abs(sp_bi[0] - c) < mach_eps) for sp_bi in sp_b):
				# if any(isclose(c, sp_bi[0]) for sp_bi in sp_b):
				# 	continue

				assert len( future[j] ) > 0, "len(future[j]) == 0"

				# Get c_dn for c in [a_min, a_max]
				if a_min <= c <= a_max or isclose( c, a_min ) or isclose( c, a_max ):
					c_dn = ((i + 1) * c + s0 * u ** (-i + 2 * j - 1)) / (i + 2)
					future_node = future[j]

					v_c_dn = 0
					if isclose( c_dn, future_node[-1][0] ):
						v_c_dn = future_node[-1][1]
					elif isclose( c_dn, future_node[0][0] ):
						v_c_dn = future_node[0][1]
					else:
						for kc in range( len( future_node ) - 1 ):
							if future_node[kc][0] <= c_dn < future_node[kc + 1][0]:
								v_c_dn = ((future_node[kc + 1][1] - future_node[kc][1]) /
								          (future_node[kc + 1][0] - future_node[kc][0]) *
								          (c_dn - future_node[kc][0]) + future_node[kc][1])
								break
					sp_c.append( (c, p_u * v_a + p_d * v_c_dn) )

			# Aggregate 'B's and 'C's.
			now_node = sorted( sp_b + sp_c, key = lambda sp_l:sp_l[0] )
			
			# Remove singular points very close to each other.
			l = 0
			while l < len( now_node ) - 1:
				if isclose( now_node[l + 1][0], now_node[l][0] ):
					if ( now_node[l][0] - s0 ) == min( now_node[l][0] - s0, now_node[l+1][0] - s0 ):
						del now_node[l + 1]
					else:
						del now_node[l]
				l += 1

			# American
			if am:
				# TODO: Remove this
				try:
					1 + now_node[-1][1]
				except IndexError:
					print( 'DEBUG: SP({i},{j}) = {spij}'.format( i = i, j = j,
					                                             spij = now_node ) )
				# print('DEBUG: SP({ip},{j}) = {spipj}'.format(ip=i+1, j=j, spij=now_node))
				if a_max - k <= now_node[-1][1]:
					pass  # Same as the European case
				elif now_node[0][1] <= a_min - k:
					now_node = [(a_min, a_min - k), (a_max, a_max - k)]
				else:
					for l in range( len(
							now_node ) - 1 ):  # TODO: Replace by binary search
						if isclose( now_node[l][0], now_node[l][1] ):
							now_node = now_node[0:l + 1]
							now_node.append( (a_max, a_max - k) )
							break
						if (now_node[l][0] - k <= now_node[l][1]) and (
									now_node[l + 1][1] <= now_node[l + 1][0] - k):
							# Find the point of intersection
							# x = ((a2 - a1) * k - (a2 * v1 - a1 * v2)) / ((A2 - v2) - (A1 - v1))
							a_int = (((now_node[l + 1][0] - now_node[l][0]) * k +
							          (now_node[l + 1][0] * now_node[l][1] -
							           now_node[l][0] * now_node[l + 1][1])) /
							         ((now_node[l + 1][0] - now_node[l][0]) -
							          (now_node[l + 1][1] - now_node[l][1])))
							now_node = now_node[0:l + 1]
							now_node.extend(
								[(a_int, a_int - k), (a_max, a_max - k)] )
							break

			# Approximations
			if h > 0.: # and (i < n - 1):
				if ub:    # Lemma 1
					l = 1
					last_point_removed = False
					while l < len( now_node ) - 1:
						if last_point_removed:
							last_point_removed = False
							l += 1
							continue
						eps = ((now_node[l + 1][1] - now_node[l - 1][1]) /
						       (now_node[l + 1][0] - now_node[l - 1][0]) *
						       (now_node[l][0] - now_node[l - 1][0]) +
						       (now_node[l - 1][1] - now_node[l][1]))
						if eps < h:  # Remove this point
							del now_node[l]
							last_point_removed = True
						else:
							l += 1

				if not ub:  # Lemma 2
					l = 1
					while l < len( now_node ) - 2:
						m1 = (now_node[l][1] - now_node[l - 1][1]) / \
						     (now_node[l][0] - now_node[l - 1][0])
						m2 = (now_node[l + 2][1] - now_node[l + 1][1]) / \
						     (now_node[l + 2][0] - now_node[l + 1][0])
						x_bar = (((m2 * now_node[l + 1][0] -
						           m1 * now_node[l - 1][0]) -
						          (now_node[l + 1][1] - now_node[l - 1][1])) /
						         (m2 - m1))
						y_bar = m1 * (x_bar - now_node[l - 1][0]) + now_node[l - 1][
							1]
						dlt = ((now_node[l + 1][1] - now_node[l][1]) /
						       (now_node[l + 1][0] - now_node[l][0]) *
						       (x_bar - now_node[l][0]) +
						       now_node[l][1] - y_bar)
						if dlt < h:  # Remove l-1 and l and add x_bar
							now_node[l] = (x_bar, y_bar)
							del now_node[l + 1]
						l += 1

			now[j] = now_node

		disp = ""
		for (j, tj) in enumerate(now):
			disp = 'N({i},{j}): {tj}\n'.format(i=i, j=j, tj=tj) + disp
		disp = '-' * 64 + '\n' + disp
		display += disp

	# print(display)

	return now[0][0][1]
