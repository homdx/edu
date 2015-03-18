class Market:
	def __init__(self, r, discrete = True):
		self.r = r
		self.discrete = discrete

class Stock:
	def __init__(self, s0, sigma):
		self.s0 = s0
		self.sigma = sigma

class EuropeanOption:
	def __init__(self, stock):
		self.stock = stock
