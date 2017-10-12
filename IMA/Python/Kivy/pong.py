import kivy
kivy.require('1.10.0')    # replace with your current Kivy version!

from kivy.app import App
from kivy.uix.widget import Widget

class PongGame(Widget):
	pass

class PongApp(App):
	def build(self):
		return PongGame()

if __name__ == "__main__":
	PongApp().run()
