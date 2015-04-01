__author__ = 'Sudip Sinha'

import math

# Computer
eps = 65536 * (7 / 3 - 4 / 3 - 1)
n = 25

# Market
r = 0.1

# Underlying
s0 = 100.
sigma = 0.2
q = 0.03

# Derivative
t = 1.
k = 90.


def pos(x: float) -> float:
	"""Returns the positive part of x"""
	return x if (x > 0) else 0


def get_underlying_tree(s0: float, sigma: float, t: float, n: int) -> list:
	"""Generate the tree of stock prices"""

	s = [[0] * (i+1) for i in range(n+1)]
	s[0][0] = s0

	dt = t / n
	u = math.exp(sigma * math.sqrt(dt))

	for i in range(1, n+1):
		for j in range(i+1):
			s[i][j] = s0 * u**(-i+2*j)

	return s


# @profile
def sp_asian_e_call(s0: float, k: float, sigma: float, q: float, t: float, n: int) -> list:
	"""Prices of an European Asian call option using the singular point method"""

	dt = t / n
	R = math.exp((r - q) * dt)
	u = math.exp(sigma * math.sqrt(dt))
	d = 1 / u
	p = (R - d) / (u - d)
	sp = [[[]] * (i + 1) for i in range(n + 1)]

	# Singular point for N(n,0)
	a = s0/(n+1) * (1 - d**(n+1)) / (1-d)
	sp[n][0] = [(a, pos(a - k))]

	# Singular point for N(n,n)
	a = s0/(n+1) * (1 - u**(n+1)) / (1-u)
	sp[n][-1] = [(a, pos(a - k))]

	# Singular points for N(n,j)    [j != 0, n]
	for j in range(1, n):
		a_min = s0 / (n+1) * ((1 - d**(n-j+1)) / (1-d) + d**(n-j-1) * (1 - u**j) / (1-u))
		a_max = s0 / (n+1) * ((1 - u**(j+1)) / (1-u) + u**(j - 1) * (1 - d**(n-j)) / (1-d))
		if k < a_min:
			sp[n][j] = [(a_min, a_min - k), (a_max, a_max - k)]
		elif k > a_max:
			sp[n][j] = [(a_min, 0), (a_max, 0)]
		else:
			sp[n][j] = [(a_min, 0), (k, 0), (a_max, a_max - k)]

	# Singular points for N(i,j)    [i != n]
	for i in range(n-1, -1, -1):
		# print('DEBUG: n = {n}, N({i},*)'.format(i=i, n=n))

		# Singular point for N(i,0)
		sp[i][0] = (p * sp[i+1][1][0][1] + (1 - p) * sp[i+1][0][0][1]) / R

		# Singular point for N(i,i)
		sp[i][i] = (p * sp[i+1][i+1][0][1] + (1 - p) * sp[i+1][i][-1][1]) / R

		# Singular point for N(i,j)    [j != 0,i]
		for j in range(1, i):
			
			a_min = s0 / (i+1) * ((1 - d**(i-j+1)) / (1-d) + d**(i-j-1) * (1 - u**j) / (1-u))
			a_max = s0 / (i+1) * ((1 - u**(j+1)) / (1-u) + u**(j-1) * (1 - d**(i-j)) / (1-d))

			if a_min > a_max:
				(a_min, a_max) = (a_max, a_min)    # Jugaad

			# 'B'
			sp_b = list()    # Initialize sp_b
			for (a, v_a) in sp[i+1][j]:
				b = ((i+2) * a - s0 * u**(-i+2*j-1)) / (i+1)

				# Get b_up for b_up in [a_min, a_max]
				if a_min - eps <= b <= a_max + eps:
					b_up = ((i+1) * b + s0 * u**(-i+2*j+1)) / (i+2)
					spi = sp[i+1][j+1]

					v_b_up = 0
					for kb in range(len(spi)):
						if kb < len(spi)-1:
							if spi[kb][0] - eps < b_up < spi[kb][0] + eps:
								v_b_up = spi[kb][1]
								break
							if spi[kb][0] - eps < b_up < spi[kb + 1][0]:
								v_b_up = ((spi[kb+1][1] - spi[kb][1]) /
								          (spi[kb+1][0] - spi[kb][0]) *
								          (b_up - spi[kb][0]) + spi[kb][1])
								break
						else:
							v_b_up = spi[kb][1]
							break
					sp_b.append((b, (p * v_b_up + (1 - p) * v_a) / R))

			#  'C'
			sp_c = list()    # Initialize sp_c
			for (a, v_a) in sp[i+1][j+1]:
				c = ((i+2) * a - s0 * u**(-i+2*j+1)) / (i+1)
				
				# Uniqueness: Verify that 'C' is not in the set of 'B'
				if any((abs(sp_bi[0] - c) < eps) for sp_bi in sp_b):
					continue
				
				# Get c_dn for c_dn in [a_min, a_max]
				if a_min - eps <= c <= a_max + eps:
					c_dn = ((i+1) * c + s0 * u**(-i+2*j-1)) / (i+2)
					spi = sp[i+1][j]

					v_c_dn = 0
					for kc in range(len(spi)):
						if kc < len(spi)-1:
							if spi[kc][0] - eps < c_dn < spi[kc][0] + eps:
								v_c_dn = spi[kc][1]
								break
							if spi[kc][0] - eps < c_dn < spi[kc + 1][0]:
								v_c_dn = ((spi[kc + 1][1] - spi[kc][1]) /
								          (spi[kc + 1][0] - spi[kc][0]) *
								          (c_dn - spi[kc][0]) + spi[kc][1])
								break
						else:
							v_c_dn = spi[kc][1]
							break
					sp_c.append((c, (p * v_a + (1 - p) * v_c_dn) / R))

			# Aggregate 'B's and 'C's
			sp_eu = sorted(sp_b + sp_c, key=lambda spi: spi[0])
			sp[i][j] = sp_eu

	return sp


def sp_asian_call(s0: float, k: float, sigma: float, q: float, t: float, n: int, american: bool=True, approx: float=0.) -> list:
	"""Prices of an American Asian call option using the singular point method"""

	dt = t / n
	R = math.exp((r - q) * dt)
	u = math.exp(sigma * math.sqrt(dt))
	d = 1 / u
	p = (R - d) / (u - d)
	sp = [[[]] * (i + 1) for i in range(n + 1)]
	if approx > 0:
		h = approx / (n * n)

	# Singular points for N(n,0)
	a = s0/(n+1) * (1 - d**(n+1)) / (1-d)
	sp[n][0] = [(a, pos(a - k))]

	# Singular points for N(n,n)
	a = s0/(n+1) * (1 - u**(n+1)) / (1-u)
	sp[n][-1] = [(a, pos(a - k))]

	# Singular points for N(n,j)
	for j in range(1, n):
		a_min = s0 / (n+1) * ((1 - d**(n-j+1)) / (1-d) + d**(n-j-1) * (1 - u**j) / (1-u))
		a_max = s0 / (n+1) * ((1 - u**(j+1)) / (1-u) + u**(j - 1) * (1 - d**(n-j)) / (1-d))
		if k < a_min:
			sp[n][j] = [(a_min, a_min - k), (a_max, a_max - k)]
		elif k > a_max:
			sp[n][j] = [(a_min, 0), (a_max, 0)]
		else:
			sp[n][j] = [(a_min, 0), (k, 0), (a_max, a_max - k)]

	# Singular points for N(i,j) i != n
	for i in range(n-1, -1, -1):
		# print('DEBUG: n = {n}, N({i},*)'.format(i=i, n=n))

		for j in range(i+1):

			a_min = s0 / (i+1) * ((1 - d**(i-j+1)) / (1-d) + d**(i-j-1) * (1 - u**j) / (1-u))
			a_max = s0 / (i+1) * ((1 - u**(j+1)) / (1-u) + u**(j-1) * (1 - d**(i-j)) / (1-d))

			if a_min > a_max:
				(a_min, a_max) = (a_max, a_min)    # Jugaad

			# 'B'
			sp_b = list()    # Initialize sp_b
			for (a, v_a) in sp[i+1][j]:
				b = ((i+2) * a - s0 * u**(-i+2*j-1)) / (i+1)

				# Get b_up for b_up in [a_min, a_max]
				if a_min - eps <= b <= a_max + eps:
					b_up = ((i+1) * b + s0 * u**(-i+2*j+1)) / (i+2)
					spi = sp[i+1][j+1]

					v_b_up = 0
					for kb in range(len(spi)):
						if kb < len(spi)-1:
							if spi[kb][0] - eps < b_up < spi[kb][0] + eps:
								v_b_up = spi[kb][1]
								break
							if spi[kb][0] - eps < b_up < spi[kb + 1][0]:
								v_b_up = ((spi[kb+1][1] - spi[kb][1]) /
								          (spi[kb+1][0] - spi[kb][0]) *
								          (b_up - spi[kb][0]) + spi[kb][1])
								break
						else:
							v_b_up = spi[kb][1]
							break
					sp_b.append((b, (p * v_b_up + (1 - p) * v_a) / R))

			#  'C'
			sp_c = list()    # Initialize sp_c
			for (a, v_a) in sp[i+1][j+1]:
				c = ((i+2) * a - s0 * u**(-i+2*j+1)) / (i+1)

				# Uniqueness: Verify that 'C' is not in the set of 'B'
				if any((abs(sp_bi[0] - c) < eps) for sp_bi in sp_b):
					continue

				# Get c_dn for c_dn in [a_min, a_max]
				if a_min - eps <= c <= a_max + eps:
					c_dn = ((i+1) * c + s0 * u**(-i+2*j-1)) / (i+2)
					spi = sp[i+1][j]

					v_c_dn = 0
					for kc in range(len(spi)):
						if kc < len(spi)-1:
							if spi[kc][0] - eps < c_dn < spi[kc][0] + eps:
								v_c_dn = spi[kc][1]
								break
							if spi[kc][0] - eps < c_dn < spi[kc + 1][0]:
								v_c_dn = ((spi[kc + 1][1] - spi[kc][1]) /
								          (spi[kc + 1][0] - spi[kc][0]) *
								          (c_dn - spi[kc][0]) + spi[kc][1])
								break
						else:
							v_c_dn = spi[kc][1]
							break
					sp_c.append((c, (p * v_a + (1 - p) * v_c_dn) / R))

			# Aggregate 'B's and 'C's
			sp_eu = sorted(sp_b + sp_c, key=lambda spi: spi[0])

			# American
			if american:
				if a_max - k <= sp_eu[-1][1] + eps:
					sp[i][j] = sp_eu
				elif sp_eu[0][1] - eps <= a_min - k:
					sp[i][j] = [(a_min, a_min - k), (a_max, a_max - k)]
				else:
					for l in range(len(sp_eu)-1):    # TODO: Replace by binary search
						if sp_eu[l][1] - eps < sp_eu[l][0] - k < sp_eu[l][1] + eps:
							sp[i][j] = sp_eu[0:l+1].append((a_max, a_max - k))
							break
						if (sp_eu[l][0] - k <= sp_eu[l][1] + eps) and (sp_eu[l+1][1] - eps <= sp_eu[l+1][0] - k):
							# Find the point of intersection
							# x = ((a2 - a1) * k - (a2 * v1 - a1 * v2)) / ((A2 - v2) - (A1 - v1))
							a_int = (((sp_eu[l+1][0] - sp_eu[l][0]) * k +
							          (sp_eu[l+1][0] * sp_eu[l][1] - sp_eu[l][0] * sp_eu[l+1][1])) /
							         ((sp_eu[l+1][0] - sp_eu[l][0]) - (sp_eu[l+1][1] - sp_eu[l][1])))
							sp_eu = sp_eu[0:l+1]
							sp_eu.extend([(a_int, a_int - k), (a_max, a_max - k)])
							sp[i][j] = sp_eu
							break
			else:
				sp[i][j] = sp_eu

	return sp


def run_test(s0: float, k: float, sigma: float, q: float, t: float, n: int, d: int):
	"""Display detailed results for a single value of n."""

	# Get the tree for the underlying and singular points
	tree = get_underlying_tree(s0=s0, sigma=sigma, t=t, n=n)
	price = sp_asian_call(s0=s0, k=k, sigma=sigma, q=q, t=t, n=n, american=True)

	# Round values to 'd' decimal places
	for i in range(n+1):
		for j in range(i+1):
			tree[i][j] = round(tree[i][j], d)
			for (l, pt) in enumerate(price[i][j]):
				price[i][j][l] = (round(price[i][j][l][0], d), round(price[i][j][l][1], d))

	# Display the tree for the underlying
	print('Tree for the underlying')
	for i in range(n + 1):
		print('{i}: {t}'.format(i=i, t=tree[i]))
	print()

	# # Display singular points
	# print('Singular Points')
	# for (i, spi) in enumerate(price):
	# 	print('i =', i)
	# 	for spij in spi:
	# 		print(spij)
	# print()

	# Display the final price of the call.
	print('Price of the call at t=0: {prc}.'.format(prc=price[0][0][0][1]))


def run_all_tests(ns: range):
	"""Display short results for a list of 'n's."""

	for i in ns:
		try:
			price = sp_asian_call(s0=s0, k=k, sigma=sigma, q=q, t=t, n=i, american=True)
			print('n={n:2}: price={prc:0.6f}'.format(n=i, prc=price[0][0][0][1]))
		except ArithmeticError:
			print('Fail (ArithmeticError): n={n:2}'.format(n=i))
		except AssertionError:
			print('Fail (AssertionError): n={n:2}'.format(n=i))


# # Display detailed results for n
# run_test(s0=s0, k=k, sigma=sigma, q=q, t=t, n=n, d=6)

# Display the results for n in [25, 50, 75, 100]
run_all_tests(range(1, 33))

# sp_asian_e_call(s0=s0, k=k, sigma=sigma, q=q, t=t, n=4)
