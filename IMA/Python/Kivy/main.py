from kivy.app import App
# from kivy.uix.button import Button
# from kivy.uix.label import Label
# from kivy.uix.scatter import Scatter
# from kivy.uix.floatlayout import FloatLayout
# from kivy.uix.textinput import TextInput
from kivy.uix.boxlayout import BoxLayout

class ScatterTextWidget(BoxLayout):
	pass

class TutorialApp(App):
	# pass
	def build(self):
		# # return Button(text = "BestFit!!!", font_size = 128)
		# bl = BoxLayout(orientation = "vertical")
		# fl = FloatLayout()
		# sc = Scatter()

		# lbl = Label(text = "Hello!!!", font_size = 64)
		# ti = TextInput(text = "Team BestFit", font_size = 64, size_hint_y = None, height = 192)
		# # ti.bind(text = ti_change)
		# ti.bind(text = lbl.setter("text"))

		# bl.add_widget(ti)
		# bl.add_widget(fl)
		# fl.add_widget(sc)
		# sc.add_widget(lbl)
		# return bl
		return ScatterTextWidget()

# def ti_change(*args):
# 	print("Text changed!")

if __name__ == "__main__":
	TutorialApp().run()
