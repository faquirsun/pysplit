from PyQt5 import QtGui, QtCore
import PyQt5.QtWidgets as qt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import cartopy.crs as ccrs

class MplCanvas(FigureCanvas):
	"""Class to represent the FigureCanvas widget"""
	def __init__(self):
		# Setup Matplotlib Figure and Axis
		self.fig = Figure()
		#self.ax = self.fig.add_subplot(111)
		self.fig.set_tight_layout(True)

		# initialisation of the canvas
		FigureCanvas.__init__(self, self.fig)
		FigureCanvas.setSizePolicy(self, qt.QSizePolicy.Expanding, qt.QSizePolicy.Expanding)

		# notify the system of updated policy
		FigureCanvas.updateGeometry(self)

	def _teleseismicMap(self, lon_cen, lat_cen):
		# Create a Cartopy axis of the calculated region and draw any coastlines
		self.ax = self.fig.add_axes([0.01, 0.01, 0.98, 0.98], projection=ccrs.AzimuthalEquidistant(central_longitude=lon_cen, central_latitude=lat_cen))
		self.ax.set_global()
		self.ax.coastlines(resolution='10m')

		return self.ax

	def _localMap(self, lon0, lat0, lon1, lat1):
		# Create a Cartopy axis of the region and draw any coastlines
		self.ax = self.fig.add_axes([0.001, 0.001, 0.998, 0.998], projection=ccrs.PlateCarree())
		#self.ax.set_global()
		self.ax.coastlines(resolution='10m')

		return self.ax

	def _tracePlot(self):
		self.ax = self.fig.add_subplot(111)


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
