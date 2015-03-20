__author__ = 'Sudip Sinha'

import math

# Computer
eps = 256 * (7 / 3 - 4 / 3 - 1)
n = 4
# TODO: Does not work for 7, 8, 16, 20, 23, 25
# TODO: Weird result for 13, 14, 17, 18, 19, 21, 22, 24

# Market
r = 0.05

# Underlying
s0 = 100.0
sigma = 0.5
q = 0.0

# Derivative
t = 0.25
k = 100.0


def pos(x):
	"""
	Returns the positive part of x
	:param x: float
	:return: float
	"""
	return x if (x > 0) else 0


def get_underlying_tree(s0, sigma, n, t):
	"""
	Generate the tree of stock prices
	:type s0: float
	:type sigma: float
	:type n: int
	:type t: float
	"""

	s = [[0] * (i + 1) for i in range(n + 1)]
	s[0][0] = s0

	dt = t / n
	u = math.exp(sigma * math.sqrt(dt))

	for i in range(1, n + 1):
		for j in range(i + 1):
			s[i][j] = s0 * u ** (-i + 2*j)

	return s


def sp_asian_e_call(s0, k, sigma, q, n, t):
	"""
	:returns Prices of an Asian call option using the singular point method.
	:rtype : list
	:type s0: float
	:type k: float
	:type sigma: float
	:type q: float
	:type n: int
	:type t: float
	"""

	dt = t / n
	u = math.exp(sigma * math.sqrt(dt))
	d = 1 / u
	R = math.exp((r - q) * dt)
	p = (R - d) / (u - d)

	# print("DEBUG: u*d = {ud}".format(ud=u*d))
	print("DEBUG: u = {u}, p = {p}, R = {R}".format(u=u, p=p, R=R))
	print()

	# Obtain the binomial tree
	sp = [[[]] * (i + 1) for i in range(n + 1)]

	# Averages for the nodes at maturity (i = n)
	# ------------------------------------------

	#  Average for N(n,0)
	a = s0/(n+1) * (1 - d**(n+1)) / (1-d)
	sp[n][0] = [(a, pos(a - k))]

	#  Average for N(n,n)
	a = s0/(n+1) * (1 - u**(n+1)) / (1-u)
	sp[n][-1] = [(a, pos(a - k))]

	# Averages for all the other nodes
	for j in range(1, n):
		a_min = s0/(n+1) * ((1 - d**(n-j+1)) / (1-d) + d**(n-j-1) * (1 - u**j) / (1-u))
		a_max = s0/(n+1) * ((1 - u**(j+1)) / (1-u) + u**(j-1) * (1 - d**(n-j)) / (1-d))
		if k < a_min:
			sp[n][j] = [(a_min, a_min - k), (a_max, a_max - k)]
		elif k > a_max:
			sp[n][j] = [(a_min, 0), (a_max, 0)]
		else:
			sp[n][j] = [(a_min, 0), (k, 0), (a_max, a_max - k)]
	# print("DEBUG: sp[maturity] = {sp}".format(sp=sp[n]))
	# print()

	# Averages for the nodes (i,j) for all i != n
	#  ------------------------------------------
	for i in range(n - 1, -1, -1):

		#  SP for all other nodes N(i,j) s.t. j != 0,i
		for j in range(0, i + 1):

			# print("DEBUG: N[{i},{j}]: #SP(N({ipp},{j}))={n_ipp_j}, #SP(N({ipp},{jpp}))={n_ipp_jpp}.".
			# 	      format(i=i, j=j, ipp=i+1, jpp=j+1, n_ipp_j=len(sp[i+1][j]), n_ipp_jpp=len(sp[i+1][j+1])))

			a_min = s0/(i+1) * ((1 - d**(i-j+1)) / (1-d) + d**(i-j-1) * (1 - u**j) / (1-u))
			a_max = s0/(i+1) * ((1 - u**(j+1)) / (1-u) + u**(j-1) * (1 - d**(i-j)) / (1-d))

			#  'B'
			# Initialize sp_b
			sp_b = list()
			for (a, v_a) in sp[i+1][j]:
				b = ((i+2) * a - s0 * u**(-i+2*j-1)) / (i+1)
				b_up = ((i+1) * b + s0 * u**(-i+2*j+1)) / (i+2)
				up_pts = sp[i+1][j+1]

				# Get b_up for b_up in [a_min, a_max]
				if a_min - eps <= b <= a_max + eps:
					v_b_up = 0
					for kb in range(len(up_pts)):
						if kb < len(up_pts)-1:
							if up_pts[kb][0] - eps < b_up < up_pts[kb][0] + eps:
								v_b_up = up_pts[kb][1]
								break
							if up_pts[kb][0] - eps < b_up < up_pts[kb+1][0]:
								v_b_up = (up_pts[kb+1][1] - up_pts[kb][1]) / \
								         (up_pts[kb+1][0] - up_pts[kb][0]) * \
								         (b_up - up_pts[kb][0]) + up_pts[kb][1]
								break
						else:
							v_b_up = up_pts[kb][1]
							break
					sp_b.append((b, 1/R * (p * v_b_up + (1 - p) * v_a)))
					# print("Inside 'B': b_up = {b_up}, up_pts = {up_pts}".format(b_up=b_up, up_pts=up_pts))

			#  'C'
			#  Initialize sp_c
			sp_c = list()
			for (a, v_a) in sp[i+1][j+1]:
				c = ((i+2) * a - s0 * u**(-i+2*j+1)) / (i+1)
				c_dn = ((i+1) * c + s0 * u**(-i+2*j-1)) / (i+2)
				dn_pts = sp[i+1][j]

				# Get c_dn for c_dn in [a_min, a_max]
				if a_min - eps <= c <= a_max + eps:
					v_c_dn = 0
					for kc in range(len(dn_pts)):
						if kc < len(dn_pts)-1:
							if dn_pts[kc][0] - eps < c_dn < dn_pts[kc][0] + eps:
								v_c_dn = dn_pts[kc][1]
								break
							if dn_pts[kc][0] - eps < c_dn < dn_pts[kc+1][0]:
								v_c_dn = (dn_pts[kc+1][1] - dn_pts[kc][1]) / \
								         (dn_pts[kc+1][0] - dn_pts[kc][0]) * \
								         (c_dn - dn_pts[kc][0]) + dn_pts[kc][1]
								break
						else:
							v_c_dn = dn_pts[kc][1]
							break
					sp_c.append((c, 1/R * (p * v_a + (1 - p) * v_c_dn)))
					# print("Inside 'C': c_dn = {c_dn}, v_c_dn = {v_c_dn}, dn_pts = {dn_pts}"
					# 	      .format(c_dn=c_dn, v_c_dn=v_c_dn, dn_pts=dn_pts))

			# Aggregate 'B's and 'C's
			sp_raw = sorted(sp_b + sp_c, key=lambda spi: spi[0])

			# Enforce uniqueness
			unique = False
			while not unique:
				unique = True
				for l in range(len(sp_raw)-1):
					if math.fabs(sp_raw[l][0] - sp_raw[l+1][0]) < eps:
						del sp_raw[l+1]
						unique = False
						break
			sp[i][j] = sp_raw

	return sp


# Get the tree for the underlying and singular points
tree = get_underlying_tree(s0=s0, sigma=sigma, n=n, t=t)
eac = sp_asian_e_call(s0=s0, k=k, sigma=sigma, q=q, n=n, t=t)

# Display the tree for the underlying
print("Tree for the underlying")
for i in range(n + 1):
	print("{i}: {t}".format(i=i, t=tree[i]))
print()

# Display singular points
print("Singular Points")
for (i, ptj) in enumerate(eac):
	# if i > 4:
	# 	continue
	print("i =", i)
	for j in ptj:
		print(j)
print()

print("Price of the call at t=0: {prc}.".format(prc=eac[0][0][0][1]))
