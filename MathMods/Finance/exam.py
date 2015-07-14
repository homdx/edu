from math import exp, log, sqrt, erf, pi

n = 2
t = 6 / 12
s0 = 50
sigma = 0.25
r = 1.0125	# This is capital R


def crr(n, t, s0, r, sigma):
	dt = t / n
	#if sigma != 0.:
		#u = r * exp( sigma * sqrt(dt))
		#d = r * exp(-sigma * sqrt(dt))
	(u, d) = (1.06, 0.95)
	p = (r - d) / (u - d)
	undlyn = ''
	eu_str = ''
	am_str = ''

	s = [[-1.] * (n+1) for i in range(2*n+1)]
	s[n][0] = s0

	eu = [[-1.] * (n+1) for i in range(2*n+1)]
	am = [[-1.] * (n+1) for i in range(2*n+1)]

	for i in range(1, n+1):
		for j in range(0, i+1):
			s[n-i+2*j][i] = s0 * d**j * u**(i-j)

	for j in range (0, n+1):
		eu[n-i+2*j][i] = (s[n-i+2*j][i] - 51) * (s[n-i+2*j][i] > 51)	# = max((s[n-i+2*j][i] - 51), 0)
		am[n-i+2*j][i] = (s[n-i+2*j][i] - 51) * (s[n-i+2*j][i] > 51)
	for i in range(n-1, -1, -1):
		for j in range (0, i+1):
			eu[n-i+2*j][i] = (p * eu[n-i+2*j-1][i+1] + (1-p) * eu[n-i+2*j+1][i+1]) / r
			am[n-i+2*j][i] = max((p * eu[n-i+2*j-1][i+1] + (1-p) * eu[n-i+2*j+1][i+1]) / r,
					     (s[n-i+2*j][i] - 51) * (s[n-i+2*j][i] > 51))

	for row in s:
		for val in row:
			if val >= 0.:
				# For the number of digits, use {5, ..., 12}
				undlyn += '{val:.2f}'.format(val=val) + '\t'
			else:
				undlyn += 1 * '\t'
		undlyn += '\n'
	for row in eu:
		for val in row:
			if val >= 0.:
				# For the number of digits, use {5, ..., 12}
				eu_str += '{val:.3f}'.format(val=val) + '\t'
			else:
				eu_str += 1 * '\t'
		eu_str += '\n'
	for row in am:
		for val in row:
			if val >= 0.:
				# For the number of digits, use {5, ..., 12}
				am_str += '{val:.3f}'.format(val=val) + '\t'
			else:
				am_str += 1 * '\t'
		am_str += '\n'

	print('Underlying\n' + undlyn)
	print('European\n' + eu_str)
	print('American\n' + am_str)
	print('u = {u:.6f}, d = {d:.6f}, p = {p:.6f}'.format(u=u, d=d, p=p))
crr(n, t, s0, r, sigma)


def N(x: float):
	return (1 + erf(x / sqrt(2))) / 2

def bs(t, x, k, r, sigma):
	d1 = log(x / k) + (r + sigma*sigma/2) * t
	d2 = log(x / k) + (r - sigma*sigma/2) * t
	pr = x * N(d1) - k * exp(-r * t) * N(d2)

	# Greeks
	delta = N(d1)
	gamma = exp(d1*d1 / 2) / (x * sigma * sqrt(2 * pi * t))
	theta = - k * exp(-r * t) * (r * N(d2) + sigma * exp(d2*d2 / 2) / (2 * sqrt(2 * pi * t)))    # c/t
	vega = k * exp(-r * t) * t * exp(d2*d2 / 2) / sqrt(2 * pi)    # c/sigma
	
	print('Call price = {pr:.6f}'.format(pr=pr))
	print('Greeks: delta={delta:.6f}, gamma={gamma:.6f}, theta={theta:.6f}'
	      .format(delta=delta, gamma=gamma, theta=theta))
bs(t, s0, 3700, r, sigma)
