from PyQt5 import QtGui, QtCore
import PyQt5.QtWidgets as qt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class MplCanvas(FigureCanvas):
	"""Class to represent the FigureCanvas widget"""
	def __init__(self):
		# Setup Matplotlib Figure and Axis
		self.fig = Figure()
		self.ax = self.fig.add_subplot(111)
		self.fig.set_tight_layout(True)

		# initialisation of the canvas
		FigureCanvas.__init__(self, self.fig)
		FigureCanvas.setSizePolicy(self, qt.QSizePolicy.Expanding, qt.QSizePolicy.Expanding)

		# notify the system of updated policy
		FigureCanvas.updateGeometry(self)


class MplWidget(qt.QWidget):
	"""Widget dedfined in Qt Designer"""
	def __init__(self, parent=None):
		# initialisation of Qt MainWindow widget
		qt.QWidget.__init__(self, parent)
		# set the canvas to the Matplotlib widget
		self.canvas = MplCanvas()
		# create a vertical box layout
		self.vbl = qt.QVBoxLayout()
		# add mpl widget to vertical box
		self.vbl.addWidget(self.canvas)
		# set the layout to the vertical box
		self.setLayout(self.vbl)
