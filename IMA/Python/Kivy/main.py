from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
# from kivy.properties import ListProperty
from random import random

EXERCISES = ["SQ", "PU", "SU"]

class DisplayWindow(BoxLayout):

	def press_playpause(self, *args):
		id_pp = self.ids.button_playpause
		if id_pp.text == "Start":
			id_pp.text = "Pause"
			# Start recording
		elif id_pp.text == "Continue":
			id_pp.text = "Pause"
			# Insert data into the table
		else:
			id_pp.text = "Continue"
			# Pause recording
		
		# Update
		exrcs = EXERCISES[int(random() * 3)]
		cnt = str(int(random() * 10))
		self.update(exercise = exrcs, count = cnt)

	def press_stop(self, *args):
		self.ids.button_playpause.text = "Start"
		
		# Update
		exrcs = EXERCISES[int(random() * 3)]
		cnt = str(int(random() * 10))
		self.update(exercise = exrcs, count = cnt)

	def update(self, exercise, count, *args):
		self.ids.exercise.text = exercise
		self.ids.count.text = count


class BestFitApp(App):
	"""docstring for BestFitApp"""
	def build(self):
		return DisplayWindow()

	
if __name__ == "__main__":
	BestFitApp().run()
