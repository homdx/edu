__author__ = 'Sudip Sinha'

import math

# Computer
eps = 65536 * (7 / 3 - 4 / 3 - 1)
n = 8

# Market
r = 0.1

# Underlying
s0 = 100.0
sigma = 0.1
q = 0.0

# Derivative
t = 0.25
k = 100.0


def pos(x: float) -> float:
	"""Returns the positive part of x"""
	return x if (x > 0) else 0


def get_underlying_tree(s0: float, sigma: float, n: int, t: float) -> list:
	"""Generate the tree of stock prices"""

	s = [[0] * (i+1) for i in range(n+1)]
	s[0][0] = s0

	dt = t / n
	u = math.exp(sigma * math.sqrt(dt))

	for i in range(1, n+1):
		for j in range(i+1):
			s[i][j] = s0 * u**(-i+2*j)

	return s


def sp_asian_e_call(s0: float, k: float, sigma: float, q: float, n: int, t: float) -> list:
	"""Prices of an Asian call option using the singular point method"""

	dt = t / n
	R = math.exp((r - q) * dt)
	u = math.exp(sigma * math.sqrt(dt))
	d = 1 / u
	p = (R - d) / (u - d)
	sp = [[[]] * (i + 1) for i in range(n + 1)]

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
			sp[i][j] = sorted(sp_b + sp_c, key=lambda spi: spi[0])

	return sp
