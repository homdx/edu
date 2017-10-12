from math import pi, ceil, log
from time import sleep

def area(width: int, height: int):
	return width * height


def n_terms_approx_pi(epsilon: float):
	n = 1
	approx = 1.
	while abs(approx * 4 - pi) > epsilon:
		n = n + 1
		approx = approx - (-1)**n / (2 * n - 1)
		# print("DEBUG: n = {n:6}, error = {err:.8f}".format(n = n, err = abs(approx * 4 - pi)))
	
	return n


# for i in range(10):
# 	print("epsilon = 10^-{i}, terms = {terms}".format(i = i, terms = n_terms_approx_pi(10**(-i))))

# epsilon = 10^-0, terms = 1
# epsilon = 10^-1, terms = 10
# epsilon = 10^-2, terms = 100
# epsilon = 10^-3, terms = 1000
# epsilon = 10^-4, terms = 10000
# epsilon = 10^-5, terms = 100001
# epsilon = 10^-6, terms = 1000001
# epsilon = 10^-7, terms = 10000001
# epsilon = 10^-8, terms = 99995329
# epsilon = 10^-9, terms = 998280591


def syr(x: int):
	seq = str(x)
	while x != 1:
		if x % 2 == 0:
			x = x // 2
		else:
			x = 3 * x + 1
		seq = seq + ", " + str(x)
	return seq

# print(syr(24))

# print("reverse"[::-1])    # Look up extended slices


def max(l: list):
	m = l[0]
	for item in list:
		if item > m:
			m = item
	return m


def binary_search(sl: list, item: int):
	low, high = 0, len(sl) - 1
	mid = (low + high) // 2
	while low <= high:
		print("({l}, {m}, {h})".format(l = low, m = mid, h = high))
		if item == sl[mid]:
			return mid
		elif item < sl[mid]:
			high = mid - 1
		else:
			low = mid + 1
		mid = (low + high) // 2
		sleep(0.1)
	return - len(sl) - 1

# print(binary_search(list(range(10)), 9))


def count_duplicates(sl1: list, sl2: list):
	n = 0
	(i1, i2) = (0, 0)
	while i1 < len(sl1) and i2 < len(sl2):
		if sl1[i1] == sl2[i2]:
			n += 1
			(i1, i2) = (i1 + 1, i2 + 1)
		elif sl1[i1] < sl2[i2]:
			i1 += 1
		else:
			i2 += 1
	return n

# print(count_duplicates([1, 2, 3], [0, 1, 3, 5]))


def num_bits(n: int):
	return ceil(log(n, 2))

# print(num_bits(9))


def down_triangle(symbol, dim) -> str:
	""" Print a triangle with the given dimensions. """
	accumulator = ""

	def down_triangle_aux(symbol: str, dim: int, accumulator: str) -> str:
		""" This is the auxillary function for tail recursion. """
	
		if dim == 0:
			return accumulator
		return down_triangle_aux(symbol, dim - 1, (symbol + " ") * dim + "\n" + accumulator)

	return down_triangle_aux(symbol, dim, accumulator)

print(down_triangle('*', 4))

# https://en.wikipedia.org/wiki/Exponentiation_by_squaring