__author__ = 'Sudip Sinha'


from tr_crr import tr_underlying
from tr_vanilla import vanilla_call
from tr_asian_singularpoints import sp_asian_call


# Market
r = 0.1

# Underlying
s0 = 100.
sigma = 0.2
q = 0.03

# Derivative
t = 1.
k = 90.
am = True

# Computation
mach_eps = 65536 * (7 / 3 - 4 / 3 - 1)
n = 2
h = 0.
ub = True


def show_results() -> None:
	tree = tr_underlying(s0=s0, sigma=sigma, t=t, n=n)
	eu = vanilla_call(r=r, s0=s0, sigma=sigma, q=q, k=k, t=t, am=am, n=n)
	sp = sp_asian_call(r=r, s0=s0, sigma=sigma, q=q, k=k, t=t, am=True, n=n)

	print('n = {n}'.format(n=n))

	print('\n' + '-' * 16 + '\tUnderlying tree \t' + '-' * 16)
	# Display the tree for the underlying.
	for (i, ti) in enumerate(tree):
		print('N({i}): {ti}'.format(i=i, ti=ti))

	print('\n' + '-' * 16 + '\tEuropean call   \t' + '-' * 16)
	# Display the tree for the European call.
	for (i, ti) in enumerate(eu):
		print('N({i}): {ti}'.format(i=i, ti=ti))

	print('\n' + '-' * 16 + '\tSingular points \t' + '-' * 16)
	# Print the price of the Asian call.
	print('Price of the call = {prc}'.format(prc=sp[0][0][0][1]))


def run(ns: list, d: int=6) -> None:
	"""Display short results for a list of 'n's."""

	for n in ns:
		args['n'] = n
		ub = sp_asian_call(**args)[0][0][0][1]
		lb = sp_asian_call(**args)[0][0][0][1]
		if n < 25:
			ac = sp_asian_call(**args)[0][0][0][1]
			print('n={n:3}: {lb:.{d}f} <= {ac:.{d}f} <= {ub:.{d}f}'.format(d=d, n=n, lb=lb, ac=ac, ub=ub))
		else:
			print('n={n:3}: {lb:.{d}f} <= {ub:.{d}f}'.format(d=d, n=n, lb=lb, ub=ub))


# Display the results
# run([2, 3, 5], d=6)

# show_results()