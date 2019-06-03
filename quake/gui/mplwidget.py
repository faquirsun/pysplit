# -*- coding: utf-8 -*-

"""
This module requires that "matplotlib" and "PyQT5" are installed in your python
environment.

Widget for plotting figures in a PyQT5 GUI.

Author: Conor Bacon
"""

from PyQt5 import QtCore
import PyQt5.QtWidgets as qt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class MplCanvas(FigureCanvas):
    """
    Class to represent the FigureCanvas widget

    Attributes
    ----------
    fig : matplotlib Figure object
        Figure that contains the Axes objects on which plots are drawn

    """

    def __init__(self):
        """
        Initialisation of Qt
        """

        # self.fig = Figure()
        self.fig = plt.figure()
        self.fig.set_facecolor("#19232D")

        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, qt.QSizePolicy.Expanding,
                                   qt.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def mapPlot(self, proj):
        self.proj = proj
        # Create a Cartopy axis of the calculated region and draw coastlines
        self.ax = self.fig.add_axes([0.01, 0.01, 0.98, 0.98], projection=proj)
        self.ax.set_global()
        self.ax.stock_img()
        self.ax.coastlines(resolution='10m')

        return self.ax

    def tracePlot(self, labels=None):

        self.fig.set_tight_layout(True)

        self.axes = self.fig.subplots(3, sharex=True)
        self.fig.subplots_adjust(hspace=0)
        # self.ax = self.fig.add_subplot(111)
        for ax in self.axes:
            ax.set_facecolor("#19232D")
            ax.xaxis.label.set_color("white")
            ax.yaxis.label.set_color("white")
            ax.tick_params(axis="x", colors="white")
            ax.tick_params(axis="y", colors="white")

            if labels is not None:
                ax.set_xlabel(labels["x"])
                ax.set_ylabel(labels["y"])

            ax.set_aspect("auto")


class MplWidget(qt.QWidget):
    """
    Widget defined in Qt Designer

    Parameters
    ----------
    canvas : MplCanvas object
        A custom MplCanvas object that contains the Axes objects on which
        plots are drawn
    vbl : qt.QVBoxLayout object
        Vertical box layout object to contain the canvas object

    """

    def __init__(self, parent=None):
        """
        Initialisation of Qt MainWindow widget

        """

        qt.QWidget.__init__(self, parent)

        self.canvas = MplCanvas()

        self.canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.canvas.setFocus()

        # create a vertical box layout
        self.vbl = qt.QVBoxLayout()
        # add mpl widget to vertical box
        self.vbl.addWidget(self.canvas)
        # set the layout to the vertical box
        self.setLayout(self.vbl)
