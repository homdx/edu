from random import getrandbits

class Board:

	def __init__(self, dim: int):
		""" Initialization. """

		self.dim = dim
		self.board = [[bool(getrandbits(1)) for j in range(dim)] for i in range(dim)]
	

	def neighbors(self, row: int, col: int) -> list:
		""" Return the neighbors of (row, col). """
		
		# Vertices
		if row == 1 and col == 1:
			return [[row, col + 1], [row + 1, col]]
		elif row == 1 and col == self.dim:
			return [[row, col - 1], [row + 1, col]]
		elif row == self.dim and col == 1:
			return [[row, col + 1], [row - 1, col]]
		elif row == self.dim and col == self.dim:
			return [[row, col - 1], [row - 1, col]]
		# Edges
		elif row == 1:
			return [[row, col - 1], [row, col + 1], [row + 1, col]]
		elif row == self.dim:
			return [[row, col - 1], [row, col + 1], [row - 1, col]]
		elif col == 1:
			return [[row, col + 1], [row - 1, col], [row + 1, col]]
		elif col == self.dim:
			return [[row, col - 1], [row - 1, col], [row + 1, col]]
		# Center~
		else:
			return [[row - 1, col], [row + 1, col], [row, col - 1], [row, col + 1]]
	

	def toggle(self, row: int, col: int):
		""" Toggle the state of (row, col) and its neighbors. """
		
		nbrs = self.neighbors(row, col)
		self.board[row - 1][col - 1] = not self.board[row - 1][col - 1]
		for nbr in nbrs:
			self.board[nbr[0] - 1][nbr[1] - 1] = not self.board[nbr[0] - 1][nbr[1] - 1]

	
	def __repr__(self) -> str:
		""" Visual representation of the state of the board. """
		
		repr = ""
		for row in self.board:
			for element in row:
				if element:
					repr = repr + "o "
				else:
					repr = repr + "@ "
			repr = repr + "\n"
		return repr


	def check(self) -> bool:
		""" Check if the lights are off. """

		return all([all(row) for row in self.board])


def play():
	""" Play the game. """

	dim = int(input("What size game do you want to start with? "))
	game = Board(dim)
	print(game)
	play = True

	while play:
		row = int(input("Enter row: "))
		col = int(input("Enter column: "))
		game.toggle(row, col)
		print(game)
		if game.check():
			print("You won!!!")
			break
		play = {"Y": True, "y": True, "N": False, "n": False}[input("Do you want to keep going (Y/N)? ")]
