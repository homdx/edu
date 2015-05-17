__author__ = 'Sudip Sinha'


from tr_crr import tr_underlying
from tr_vanilla import vanilla_call
from tr_asian_singularpoints import sp_asian_call


# Market
r = 0.05

# Underlying
s0 = 50.
sigma = 0.2
q = 0.0

# Derivative
t = 1.
k = 60.
am = False

# Computation
mach_eps = 65536 * (7 / 3 - 4 / 3 - 1)
n = 5
h = 1e-4
ub = True


def show_results(n: int) -> None:
	tree = tr_underlying(s0=s0, sigma=sigma, t=t, n=n)
	eu = vanilla_call(r=r, s0=s0, sigma=sigma, q=q, k=k, t=t, am=am, n=n)
	# sp = sp_asian_call(r=r, s0=s0, sigma=sigma, q=q, k=k, t=t, am=True, n=n)

	print('n = {n}'.format(n=n))

	print('\n' + '-' * 16 + '\tUnderlying tree \t' + '-' * 16)
	# Display the tree for the underlying.
	for (i, ti) in enumerate(tree):
		print('N({i}): {ti}'.format(i=i, ti=ti))

	print('\n' + '-' * 16 + '\tEuropean call   \t' + '-' * 16)
	# Display the tree for the European call.
	for (i, ti) in enumerate(eu):
		print('N({i}): {ti}'.format(i=i, ti=ti))

	# print('\n' + '-' * 16 + '\tSingular points \t' + '-' * 16)
	# # Print the price of the Asian call.
	# for (i, ti) in enumerate(sp):
	# 	for (j, tj) in enumerate(ti):
	# 		print('N({i},{j}): {tj}'.format(i=i, j=j, tj=tj))
	# 	print('-' * 64)
	#
	# print('Price of the call = {prc}'.format(prc=sp[0][0][0][1]))


def run_asian(ns: list, d: int=6) -> None:
	"""Display short results for a list of 'n's."""

	for n in ns:
		upper = sp_asian_call(r=r, s0=s0, sigma=sigma, q=q, k=k, t=t, am=True, n=n, h=h, ub=True)[0][0][0][1]
		lower = sp_asian_call(r=r, s0=s0, sigma=sigma, q=q, k=k, t=t, am=True, n=n, h=h, ub=False)[0][0][0][1]
		if n < 25:
			actual = sp_asian_call(r=r, s0=s0, sigma=sigma, q=q, k=k, t=t, am=True, n=n)[0][0][0][1]
			print('n={n:3}: {lb:.{d}f} <= {ac:.{d}f} <= {ub:.{d}f}'.format(d=d, n=n, lb=lower, ac=actual, ub=upper))
		else:
			print('n={n:3}: {lb:.{d}f} <= {ub:.{d}f}'.format(d=d, n=n, lb=lower, ub=upper))


# Display the results
# run_asian([100], d=6)
# run_asian([221], d=6)
# run_asian([222], d=6)

show_results(n)

# print(sp_asian_call(r=r, s0=s0, sigma=sigma, q=q, k=k, t=t, am=True, n=222, h=h, ub=True)[0][0][0][1])
# print(sp_asian_call(r=r, s0=s0, sigma=sigma, q=q, k=k, t=t, am=True, n=222, h=h, ub=False)[0][0][0][1])

# http://www.goddardconsulting.ca/matlab-binomial-crr.html
# http://www.hoadley.net/options/binomialtree.aspx?tree=B
# http://www.mngt.waikato.ac.nz/kurt/frontpage/studentwork/danielchainov2003/bitree2.htm
# http://www.fintools.com/resources/online-calculators/options-calcs/binomial/
# http://www.optionspricevaluation.com/
