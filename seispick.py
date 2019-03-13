#!/usr/bin/env python3

"""
This module requires that 'numpy', 'pandas', 'matplotlib' and 'PyQT5' are
installed in your Python environment

The back-end is found within the seispick.metainfo and seispick.catalogue modules.

Author: Hemmelig
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PyQt5.QtCore import Qt, QDate, QDateTime
import PyQt5.QtWidgets as qt
import PyQt5.QtGui as QtGui
from PyQt5 import uic
import sys
import glob
import os
import re
import pathlib

import metainfo as psm
import catalogue as cat


class SeisPick(qt.QMainWindow):
	"""
	SeisPick class

	Base class that connects the PyQt5 SeisPick GUI to the back-end system of 
	classes

	Attributes
	----------

	Methods
	-------

	Properties
	----------
	"""

	def __init__(self):
		"""
		Class initialisation

		"""
		super().__init__()

		self.filt = None

		self.initUI()

	def initUI(self):
		"""
		Load UI from designer file and initialise interface

		"""

		uic.loadUi('ui_files/main.ui', self)

		self.connect()

		self.uiMapMpl.hide()
		self.uiStatusStacked.setCurrentIndex(1)

		self.setWindowTitle('SeisPick - Main catalogue window')
		self.setWindowIcon(QtGui.QIcon("misc/icon.png"))
		self.show()

	def connect(self):
		"""
		Connect all signals from Qt to functions

		"""

		self.uiNewCatalogueAction.triggered.connect(self.newCatalogue)
		self.uiLoadCatalogueAction.triggered.connect(self.loadCatalogue)
		self.uiDefaultFilterAction.triggered.connect(self.defaultFilter)

		self.uiMapMpl.canvas.mpl_connect('pick_event', self._onPick)
		self.uiMapMpl.canvas.mpl_connect('button_press_event', self._onMapClick)
		self.uiMapMpl.canvas.mpl_connect('motion_notify_event', self._onMapMove)
		self.uiMapMpl.canvas.mpl_connect('button_release_event', self._onMapRelease)

		self.uiReceiversListView.doubleClicked.connect(self.receiverSelect)
		self.uiSourcesListView.doubleClicked.connect(self.sourceSelect)
		self.uiArrivalLoadButton.clicked.connect(self.loadArrivals)
		self.uiWaveformLoadButton.clicked.connect(self.loadWaveforms)
		self.uiPlotMapButton.clicked.connect(lambda: self.plotCatalogueMap(replot=True))
		self.uiResetMapButton.clicked.connect(self.plotCatalogueMap)

	def defaultFilter(self):
		self.defaultFilterDialogue = DefaultFilterDialogue(self)

		if not self.defaultFilterDialogue.exec_():
			return

	def newCatalogue(self):
		"""
		Opens a NewCatalogueDialogue window and creates a new catalogue.

		If the dialogue is not successfully	executed, no catalogue is generated. 
		Testing of correct input is handled	in the NewCatalogueDialogue class.

		"""

		newCatalogueDialogue = NewCatalogueDialogue()

		if not newCatalogueDialogue.exec_():
			return

		catalogue_parameters = newCatalogueDialogue.catalogue_parameters

		self.createCatalogue(params=catalogue_parameters, new=True)

	def loadCatalogue(self):
		"""
		Opens a QFileDialog window and loads an existing catalogue.

		A test is established to ensure the selected directory is compatible
		with SeisPick.

		"""

		catalogue_path = qt.QFileDialog.getExistingDirectory(self, 'Choose catalogue directory')

		data_dir = pathlib.Path(catalogue_path) / "data"
		meta_dir = pathlib.Path(catalogue_path) / "metafiles"

		if not data_dir.exists() and meta_dir.exists():
			qt.QMessageBox.about(self, "Error!", "The directory you have chosen is invalid, please try again.")
			return

		self.uiStatusBar.showMessage('Parsing catalogue metafile...')
		try:
			catalogue_parameters = cat.parseCatalogueMetafile(meta_dir)
		except FileNotFoundError:
			qt.QMessageBox.about(self, "Error!", "Unable to find file containing the metadata for this catalogue. Please double-check it exists.")
			return

		self.createCatalogue(params=catalogue_parameters)

	def createCatalogue(self, params=None, new=False):
		"""
		Creates a Catalogue object from the catalogue parameter dictionary.

		Parameters
		----------
		params : dict
			Dictionary object that contains the parameters that define the 
			catalogue
		new : bool, optional
			If True, a new catalogue will be generated, including generation
			of all relevant directories and files (default: False)

		TO-DO
		-----
		If a new catalogue instance is opened in the same SeisPick session,
		need to create a function that will restore the interface to the right
		state and reset all internal variables to default values.

		"""

		self.uiCatalogueTypeDisplay.setText(params["catalogue_type"])
		self.uiCatalogueNameDisplay.setText(params["catalogue_name"])
		self.uiStatusStacked.setCurrentIndex(0)

		self.uiReceiversListView.setEnabled(True)
		self.uiSourcesListView.setEnabled(True)

		if params["catalogue_type"] == "teleseismic":
			self.uiStatusBar.showMessage("Creating an instance of TeleseismicCatalogue...")
			self.uiPlotOptionsStacked.setCurrentIndex(1)

			self.catalogue = cat.TeleseismicCatalogue(new, **params)

			self._populateCatalogueInformation(params)

		elif params["catalogue_type"] == "local":
			self.uiStatusBar.showMessage("Creating an instance of LocalCatalogue...")
			self.uiPlotOptionsStacked.setCurrentIndex(0)

			self.catalogue = cat.LocalCatalogue(new, **params)

		self._populateReceiverList()

		self.plotCatalogueMap()

	def loadArrivals(self):
		"""
		Loads an existing file or generates one

		"""

		arrivals = self.catalogue.loadArrivals()
		if arrivals:
			self.uiArrivalLoadButton.setEnabled(False)
			self.uiWaveformLoadButton.setEnabled(True)
		else:
			self.uiStatusBar.showMessage("The arrivals file does not exist - attempting to generate one...")

			if self.catalogue.catalogue_type == "local":
				self.catalogue.getArrivals()
			elif self.catalogue.catalogue_type == "teleseismic":
				self.input_phases = []
				
				telePhaseDialogue = TelePhaseDialogue()

				if not telePhaseDialogue.exec_():
					return

				self.input_phases = telePhaseDialogue.phases

				self.catalogue.getArrivals(phases=self.input_phases)
				self._populateReceiverList()

			self.loadArrivals()

	def loadWaveforms(self):
		"""
		Tests if the all waveform data has been loaded according to the 
		arrivals DataFrame

		"""

		waveforms = self.catalogue.loadWaveforms()
		if not waveforms:
			self.catalogue.getWaveforms()

		self.uiWaveformLoadButton.setEnabled(False)

	def _populateReceiverList(self):
		"""
		Populates the ReceiverList object with receiver names

		"""

		model = QtGui.QStandardItemModel(self.uiReceiversListView)
		for receiver in self.catalogue.network.receivers.name.values:
			item = QtGui.QStandardItem(receiver)
			model.appendRow(item)
		self.uiReceiversListView.setModel(model)

	def _populateCatalogueInformation(self, params):
		"""
		Populates the Catalogue Info widget with catalogue information

		Parameters
		----------
		params : dict
			Dictionary object that contains the parameters that define the 
			catalogue			

		"""

		self.uiStartDateDisplay.setText(params["start_date"])
		self.uiEndDateDisplay.setText(params["end_date"])
		self.uiCenLonDisplay.setText(params["clon"])
		self.uiCenLatDisplay.setText(params["clat"])
		self.uiMinRadDisplay.setText(params["minrad"])
		self.uiMaxRadDisplay.setText(params["maxrad"])
		self.uiMinMagDisplay.setText(params["minmag"])
		self.uiSourceCountDisplay.setText(self.catalogue.source_count)

	def plotCatalogueMap(self, replot=False):
		self.uiStatusBar.showMessage("Plotting catalogue map...")

		self.map_bg = None
		self.map_drag = False
		self.map_click = None

		self.uiMapMpl.show()

		if replot:
			self.uiMapMpl.canvas.ax.clear()
			lims = {"lon0": float(self.uiMinLonInput.text()),
					"lon1": float(self.uiMaxLonInput.text()),
					"lat0": float(self.uiMinLatInput.text()),
					"lat1": float(self.uiMaxLatInput.text())}

			if self.catalogue.catalogue_type == "local":
				self.catalogue.plotGeographic(map_widget=self.uiMapMpl, lims=lims, receivers=True)
			else:
				self.catalogue.plotGeographic(map_widget=self.uiMapMpl)

		else:
			if self.catalogue.catalogue_type == "local":
				self.catalogue.plotGeographic(map_widget=self.uiMapMpl, receivers=True)
				self.uiMinLonInput.setText(str(f"{self.catalogue.lon0:.5f}"))
				self.uiMaxLonInput.setText(str(f"{self.catalogue.lon1:.5f}"))
				self.uiMinLatInput.setText(str(f"{self.catalogue.lat0:.5f}"))
				self.uiMaxLatInput.setText(str(f"{self.catalogue.lat1:.5f}"))
			else:
				self.catalogue.plotGeographic(map_widget=self.uiMapMpl)

		# Connect to the map to grab background once Qt has done resizing
		self.uiMapMpl.canvas.mpl_connect('draw_event', self._mapDrawEvent)

		# Draw the canvas once Qt has done all resizing
		self.uiMapMpl.canvas.draw_idle()

		self.uiStatusBar.showMessage("Catalogue map plot complete.")

	def _mapDrawEvent(self, event):
		"""
		Grab the map background when it is drawn

		Parameters
		----------
		event : DrawEvent
			Contains information about the draw event

		"""

		self.map_bg = event.canvas.copy_from_bbox(event.canvas.ax.bbox)

	def _onPick(self, event):
		"""	
		Process picks on data in plot

		Parameters
		----------
		event : PickEvent
			Contains information about the pick event

		"""

		if event.mouseevent.key == "shift":
			return

		artist = event.artist

		xy = artist.get_offsets()
		ind = event.ind
		lon, lat = xy[ind[0]][0], xy[ind[0]][1]

		label = artist.get_label()
		if "REC" in label:
			rec = label.split(": ")[1]
			self.receiver = self.catalogue.network.lookupReceiver(rec)
			self._updateReceiverInformation(self.receiver)

			if event.mouseevent.key == "control":
				self.pickWindow = PickingWindow(self.catalogue, receiver=self.receiver, filt=self.filt)

		if "SOURCE" in label:
			src = label.split(": ")[1]
			self.source = self.catalogue.lookupSource(src)
			self._updateSourceInformation(self.source)

			if event.mouseevent.key == "control":
				self.pickWindow = PickingWindow(self.catalogue, source=self.source, filt=self.filt)

	def _onMapClick(self, event):
		"""

		"""

		if event.key != "shift":
			return
		if self.catalogue.catalogue_type == "teleseismic":
			return
		if self.map_drag:
			return
		if event.inaxes != self.uiMapMpl.canvas.ax:
			return

		# Grab the x and y positions of the click point
		xpress, ypress = event.xdata, event.ydata
		self.map_drag = True

		# Store the click point 
		self.map_click = xpress, ypress

	def _onMapMove(self, event):
		"""

		"""

		if event.key != "shift":
			return
		if self.catalogue.catalogue_type == "teleseismic":
			return
		if not self.map_drag:
			return
		if event.inaxes != self.uiMapMpl.canvas.ax:
			return

		# Grab the x and y data positions of the click point
		xpress, ypress = self.map_click

		# Grab the x and y data positions of the motion point
		xmove, ymove = event.xdata, event.ydata
		dx = xmove - xpress
		dy = ymove - ypress

		# Draw background from pixel buffer
		self.uiMapMpl.canvas.restore_region(self.map_bg)

		# Set rectangle values
		map_rectangle = Rectangle((xpress, ypress), dx, dy,
								  edgecolor='red', fill=False, transform=self.uiMapMpl.canvas.proj)
		self.uiMapMpl.canvas.ax.add_patch(map_rectangle)
		self.uiMapMpl.canvas.ax.draw_artist(map_rectangle)

		# Blit the redrawn area
		self.uiMapMpl.canvas.blit(self.uiMapMpl.canvas.ax.bbox)

	def _onMapRelease(self, event):
		"""

		"""

		if event.key != "shift":
			return
		if self.catalogue.catalogue_type == "teleseismic":
			return
		if not self.map_drag:
			return

		# Grab the x and y data positions of the click point
		xpress, ypress = self.map_click

		# Reset map click and lock variables
		self.map_click = None
		self.map_drag  = False

		# Grab the x and y data positions of the release point
		xrelease, yrelease = event.xdata, event.ydata
		dx = xrelease - xpress
		dy = yrelease - ypress

		# Draw background from pixel buffer
		self.uiMapMpl.canvas.restore_region(self.map_bg)

		# Draw the final rectangle
		map_rectangle = Rectangle((xpress, ypress), dx, dy,
								  edgecolor='red', fill=False)
		self.uiMapMpl.canvas.ax.add_patch(map_rectangle)
		self.uiMapMpl.canvas.ax.draw_artist(map_rectangle)

		# Blit the redrawn area
		self.uiMapMpl.canvas.blit(self.uiMapMpl.canvas.ax.bbox)

		# Set the text values of the lon/lat input boxes
		self.uiMinLonInput.setText(str(f"{min(xpress, xrelease):.2f}"))
		self.uiMaxLonInput.setText(str(f"{max(xpress, xrelease):.2f}"))
		self.uiMinLatInput.setText(str(f"{min(ypress, yrelease):.2f}"))
		self.uiMaxLatInput.setText(str(f"{max(ypress, yrelease):.2f}"))

	def receiverSelect(self, index):
		"""
		Updates the available sources for a given receiver

		Parameters
		----------
		index : int
			Index of the item picked from the Receivers List widget

		"""

		receiver = self.uiReceiversListView.model().data(index)

		self.receiver = self.catalogue.network.lookupReceiver(receiver)
		self.list_receiver = self.receiver

		sources = list(self.catalogue.arr_df.query('receiverid == @receiver').sourceid.values)

		data_dir = pathlib.Path(self.catalogue.catalogue_path) / "data" / receiver.upper()
		avail_sources = list(data_dir.glob("*.z"))
		avail_sources = [x.parts[-1] for x in avail_sources]

		model = QtGui.QStandardItemModel(self.uiSourcesListView)
		for source in sources:
			item = QtGui.QStandardItem("Source {}".format(source))
			if not any("source.{}".format(source) in s for s in avail_sources):
				item.setForeground(QtGui.QBrush(Qt.gray))
			model.appendRow(item)
		self.uiSourcesListView.setModel(model)

		self._updateReceiverInformation(self.receiver)

	def sourceSelect(self, index):
		"""
		Opens a picking window for a given source/receiver selection

		If the waveform files have not been downloaded, it attempts to retrieve
		them from the archive. If the file does not exist, the arrival is 
		removed from the catalogue's arrivals DataFrame.

		Parameters
		----------
		index : int
			Index of the item picked from the Sources List widget

		"""

		source = self.uiSourcesListView.model().data(index).split(" ")[1]

		self.source = self.catalogue.lookupSource(source)

		self._updateSourceInformation(self.source)

		waveform = self.catalogue.getWaveform(self.source, self.list_receiver)
		if not waveform:
			self.uiSourcesListView.model().removeRow(index.row())
			qt.QMessageBox.about(self, "Warning!", "No archive data available for this arrival.")
		else:
			self.pickWindow = PickingWindow(self.catalogue, 
											receiver=self.list_receiver, 
											source=self.source, 
											filt=self.filt)

	def _updateReceiverInformation(self, receiver):
		"""
		Populates the Receiver Info widget with receiver information

		"""

		self.uiReceiverNameDisplay.setText(receiver.station)
		self.uiReceiverLonDisplay.setText(f"{receiver.longitude:.4f}")
		self.uiReceiverLatDisplay.setText(f"{receiver.latitude:.4f}")
		self.uiReceiverElevDisplay.setText(f"{receiver.elevation:.4f}")
		self.uiReceiverDepDisplay.setText(receiver.deployment.isoformat().split("T")[0])
		self.uiReceiverRetDisplay.setText(receiver.retrieval.isoformat().split("T")[0])

	def _updateSourceInformation(self, source):
		"""
		Populates the Source Info widget with source information

		"""

		self.uiOriginDateDisplay.setText(source.otime.isoformat().split("T")[0])
		self.uiOriginTimeDisplay.setText(source.otime.isoformat().split("T")[1])
		self.uiSourceLonDisplay.setText(f"{source.longitude:.4f}")
		self.uiSourceLatDisplay.setText(f"{source.latitude:.4f}")
		self.uiSourceDepthDisplay.setText(f"{source.depth:.4f}")
		if (type(source.magnitude) == float) or (type(source.magnitude) == np.float64):
			self.uiSourceMagDisplay.setText(f"{source.magnitude:.2f}")
		else:
			self.uiSourceMagDisplay.setText(str(source.magnitude))
		self.uiSourceIDDisplay.setText(str(source.sourceid))


class NewCatalogueDialogue(qt.QDialog):

	# ------------------------------
	# Class initialisation functions
	# ------------------------------

	def __init__(self):
		super().__init__()

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/new_catalogue_dialogue.ui', self)

		# Initialise with teleseismic catalogue selected
		self.uiTeleseismicRadio.setChecked(True)
		self.uiInputTypeStacked.setCurrentIndex(0)

		# Set the default input options for local catalogues
		self.uiLocalOptionsStacked.setCurrentIndex(1)

		# Connect to actions and interactive widgets
		self.connect()

		self.setWindowTitle('SeisPick - New catalogue')
		self.setWindowIcon(QtGui.QIcon("misc/icon.png"))
		self.show()

	# ------------------------------

	# ------------------------------
	# Connection and event functions
	# ------------------------------

	def connect(self):
		# Connect to directory browsing buttons
		self.uiCatalogueDirectoryButton.clicked.connect(self.browseCatalogue)
		self.uiArchiveDirectoryButton.clicked.connect(self.browseArchive)
		self.uiReceiverFileButton.clicked.connect(self.browseReceivers)

		# Handle toggling of local/teleseismic radio buttons
		self.uiTeleseismicRadio.toggled.connect(self.teleRadioClicked)
		self.uiLocalRadio.toggled.connect(self.localRadioClicked)

		# Connect to accept/reject actions
		self.uiButtonBox.accepted.connect(self.actionAccept)
		self.uiButtonBox.rejected.connect(self.actionReject)

		# Connect to local input actions
		self.uiLocalFileButton.clicked.connect(self.browseLocalFile)
		self.uiLocalPathButton.clicked.connect(self.browseLocalPath)
		self.uiLocalInputFormatComboBox.currentIndexChanged.connect(self.inputSelect)		

	def browseCatalogue(self):
		catalogue_path = qt.QFileDialog.getExistingDirectory(self, 'Choose catalogue directory')
		self.uiCatalogueDirectoryInput.setText(catalogue_path)

	def browseArchive(self):
		archive_path = qt.QFileDialog.getExistingDirectory(self, 'Choose archive directory')
		self.uiArchiveDirectoryInput.setText(archive_path)

	def browseReceivers(self):
		filename = qt.QFileDialog.getOpenFileName(self, 'Open file')
		self.uiReceiverFileInput.setText(filename[0])

	def browseLocalFile(self):
		filename = qt.QFileDialog.getOpenFileName(self, 'Open file')
		self.uiLocalFileInput.setText(filename[0])

	def browseLocalPath(self):
		pathname = qt.QFileDialog.getExistingDirectory(self, 'Choose SeisLoc directory')
		self.uiLocalPathInput.setText(pathname)

	# ------------------------------

	# --------------------------
	# Handlers for input options
	# --------------------------

	def teleRadioClicked(self, enabled):
		if enabled:
			self.uiInputTypeStacked.setCurrentIndex(0)

	def localRadioClicked(self, enabled):
		if enabled:
			self.uiInputTypeStacked.setCurrentIndex(1)

	def inputSelect(self):
		# Handle each different type of input
		if self.uiLocalInputFormatComboBox.currentText() == "SeisLoc":
			self.uiLocalOptionsStacked.setCurrentIndex(1)
		else:
			self.uiLocalOptionsStacked.setCurrentIndex(0)

	# --------------------------

	# ------------------------------
	# Accept/Reject action overrides
	# ------------------------------

	def actionAccept(self):
		# If the inputs are accepted, set all of the parameters within the parent class
		# Test if any of the input fields have not been completed
		if (self.uiCatalogueDirectoryInput.text() == "") or (self.uiCatalogueNameInput.text() == "") or (self.uiReceiverFileInput.text() == ""):
			qt.QMessageBox.about(self, "Error!", "Please complete the form and try again.")
			return

		# Catalogue path
		catalogue_path = self.uiCatalogueDirectoryInput.text()
		if not os.path.exists(catalogue_path):
			qt.QMessageBox.about(self, "Error!", "You must provide a valid catalogue path.")
			return

		# Catalogue name
		catalogue_name = self.uiCatalogueNameInput.text()
		if not type(catalogue_name) == str:
			qt.QMessageBox.about(self, "Error!", "The catalogue name must be a String.")
			return

		# Archive path
		archive_path = self.uiArchiveDirectoryInput.text()
		if archive_path == "":
			pass
		elif not os.path.exists(archive_path):
			qt.QMessageBox.about(self, "Error!", "You must provide a valid archive path.")
			return

		archive_type = self.uiArchiveFormatComboBox.currentText()
		# SeisComp3 
		if archive_type == "SeisComp3":
			archive_format = "{year}/*/{receiver}/{comp}/*.{receiver}..{comp}.D.{year}.{jday}"
		# Cambridge Volcano Seismology archive format
		elif archive_type == "Cambridge Volcano Seismology":
			archive_format = "{year}/{jday}/*_{receiver}_{comp}"
		
		# Receiver file
		rec_file  = self.uiReceiverFileInput.text()
		if not os.path.exists(rec_file):
			qt.QMessageBox.about(self, "Error!", "You must provide a valid receiver file.")
			return

		if self.uiTeleseismicRadio.isChecked():
			# Set catalogue type and data source
			catalogue_type = "teleseismic"
			data_source    = "IRIS"

			# Create dictionary of variables specific to teleseismic catalogues
			d_specific = {'minmag': self.uiMinMagInput.text(),
						  'clon': self.uiCenLonInput.text(),
						  'clat': self.uiCenLatInput.text(),
						  'minrad': self.uiMinRadInput.text(),
				          'maxrad': self.uiMaxRadInput.text()}

		elif self.uiLocalRadio.isChecked():
			# Set catalogue type and data source
			catalogue_type = "local"
			data_source = self.uiLocalInputFormatComboBox.currentText()

			if data_source == "SeisLoc":
				local_input = self.uiLocalPathInput.text()
				if not os.path.exists(local_input):
					qt.QMessageBox.about(self, "Error!", "You must provide a valid input directory.")
					return
			else:
				local_input = self.uiLocalFileInput.text()
				if not os.path.isfile(local_input):
					qt.QMessageBox.about(self, "Error!", "You must provide a valid input file.")
					return

			# Create dictionary of variables specific to local catalogues
			d_specific = {'local_input': local_input}

		# Create catalogue parameter dictionary and merge with catalogue type specific parameters
		catalogue_parameters = {'catalogue_name': catalogue_name,
								'catalogue_type': catalogue_type,
								'catalogue_path': "{}/{}".format(catalogue_path, catalogue_name),
								'data_source': data_source,
								'creation_date': QDateTime.currentDateTime().toUTC().toString(Qt.ISODate),
								'archive_path': archive_path,
								'archive_type': archive_type,
								'archive_format': archive_format,
								'rec_file': rec_file,
								'start_date': self.uiStartDateInput.date().toString(Qt.ISODate),
								'end_date': self.uiEndDateInput.date().toString(Qt.ISODate)}
		self.catalogue_parameters = {**catalogue_parameters, **d_specific}

		# Send accept signal to Dialog
		self.accept()

	def actionReject(self):
		# Send reject signal to Dialog
		self.reject()

	# ------------------------------


class TelePhaseDialogue(qt.QDialog):

	# ------------------------------
	# Class initialisation functions
	# ------------------------------

	def __init__(self):
		super().__init__()

		self.phases = []

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/tele_phase_dialogue.ui', self)

		# Connect to actions and interactive widgets
		self.connect()

		# Populate the receiver list
		self.model = QtGui.QStandardItemModel(self.uiPhaseListView)

		self.setWindowTitle('SeisPick - Teleseismic phase selection')
		self.setWindowIcon(QtGui.QIcon("misc/icon.png"))
		self.show()

	# ------------------------------

	# ------------------------------
	# Connection and event functions
	# ------------------------------

	def connect(self):
		# Connect to the add phase button
		self.uiAddPhaseButton.clicked.connect(self.addPhase)

		# Connect to accept/reject actions
		self.uiButtonBox.accepted.connect(self.actionAccept)
		self.uiButtonBox.rejected.connect(self.actionReject)		

	# ------------------------------

	# --------------------------
	# Handlers for input options
	# --------------------------

	def addPhase(self):
		# Grab the current phase to be added
		phase = self.uiInputPhaseComboBox.currentText()

		# Check if phase already in list
		if phase in self.phases:
			pass
		# If not, add it
		else:
			self.phases.append(phase)

			# Add it to the ListView object
			item = QtGui.QStandardItem(phase)
			self.model.appendRow(item)

			# Update model
			self.uiPhaseListView.setModel(self.model)

	# --------------------------

	# ------------------------------
	# Accept/Reject action overrides
	# ------------------------------

	def actionAccept(self):
		if not self.phases:
			qt.QMessageBox.about(self, "Error!", "No phases selected.")
			return

		# Update the input phases in the parent class
		self.parent.input_phases = self.phases

		# Send accept signal to Dialog
		self.accept()

	def actionReject(self):
		# Send reject signal to Dialog
		self.reject()

	# ------------------------------	


class DefaultFilterDialogue(qt.QDialog):

	# ------------------------------
	# Class initialisation functions
	# ------------------------------

	def __init__(self, parent):
		super().__init__()

		self.parent = parent

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/default_filter_dialogue.ui', self)

		# Connect to actions and interactive widgets
		self.connect()

		self.setWindowTitle('SeisPick - Set default filter')
		self.setWindowIcon(QtGui.QIcon("misc/icon.png"))
		self.show()

	# ------------------------------

	# ------------------------------
	# Connection and event functions
	# ------------------------------

	def connect(self):
		# Connect to accept/reject actions
		self.uiButtonBox.accepted.connect(self.actionAccept)
		self.uiButtonBox.rejected.connect(self.actionReject)

	# ------------------------------

	# ------------------------------
	# Accept/Reject action overrides
	# ------------------------------

	def actionAccept(self):
		# Builds a dictionary of keyword arguments for the filter
		# function in the Source class, based on filter type
		filter_type = self.uiFilterTypeComboBox.currentText()
		corners     = int(self.uiPoleCountComboBox.currentText())
		zerophase   = self.uiZeroPhaseCheckBox.isChecked()

		filter_options = {"type": filter_type,
						  "corners": corners,
						  "zerophase": zerophase}

		if filter_type == "bandpass":
			try:
				freqmin = float(self.uiLowFreqInput.text())
			except ValueError:
				qt.QMessageBox.about(self, "Error!", "Invalid lowpass frequency.")
				return

			try:
				freqmax = float(self.uiHighFreqInput.text())
			except ValueError:
				qt.QMessageBox.about(self, "Error!", "Invalid highpass frequency.")
				return

			# Create filter options dictionary
			filter_freqs = {'freqmin': freqmin,
							  'freqmax': freqmax}

		elif filter_type == "lowpass":
			try:
				freq = float(self.uiLowFreqInput.text())
			except ValueError:
				qt.QMessageBox.about(self, "Error!", "Invalid lowpass frequency.")
				return

			# Create filter options dictionary
			filter_freqs = {'freq': freq}

		elif filter_type == "highpass":
			try:
				freq = float(self.uiHighFreqInput.text())
			except ValueError:
				qt.QMessageBox.about(self, "Error!", "Invalid highpass frequency.")
				return

			# Create filter options dictionary
			filter_freqs = {'freq': freq}

		# Update the default filter in parent class
		self.parent.filt = {**filter_options, **filter_freqs}

		# Send accept signal to Dialog
		self.accept()

	def actionReject(self):
		# Send reject signal to Dialog
		self.reject()

	# ------------------------------


class CustomPickDialogue(qt.QDialog):

	def __init__(self, parent):
		super().__init__()

		self.parent = parent

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/custom_phase_dialogue.ui', self)

		self.setWindowTitle('SeisPick - Set custom phase pick')
		self.setWindowIcon(QtGui.QIcon("misc/icon.png"))
		self.show()

	def actionAccept(self):
		# Read the custom pick type
		try:
			self.parent.pick_type = self.customPickDialogue.uiPhaseComboBox.currentText()
			self.uiCustomPhaseDisplay.setText(self.parent.pick_type)

			self.pick_line_color = "orange"
		except ValueError:
			qt.QMessageBox.about(self, "Error!", "You need to specify a phase to pick!")

		# Send accept signal to Dialog
		self.accept()

	def actionReject(self):
		"""
		
		"""

		self.reject()


class WadatiWindow(qt.QMainWindow):

	wad_labels  = {"x": "P traveltime / s",
			       "y": "S - P traveltime / s"}

	dist_labels = {"x": "Epicentral distance / km",
				   "y": "Traveltime / s"}

	def __init__(self, parent):
		super().__init__()

		self.ptravels  = []
		self.stravels  = []
		self.sptimes   = []
		self.receivers = []
		self.dists     = []

		self.parent = parent

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/wadati_window.ui', self)

		self.uiWadatiPlotMpl.canvas.tracePlot(labels=self.wad_labels)
		self.uiDistancePlotMpl.canvas.tracePlot(labels=self.dist_labels)

		self.connect()

		self.setWindowTitle('SeisPick - Wadati plot window')
		self.setWindowIcon(QtGui.QIcon("misc/icon.png"))
		self.show()

	def connect(self):
		# Calculate Vp/Vs button
		self.uiCalcVpVsButton.clicked.connect(self.calcVpVs)

		# Plot connection
		self.uiWadatiPlotMpl.canvas.mpl_connect('pick_event', self._onPick)
		self.uiDistancePlotMpl.canvas.mpl_connect('pick_event', self._onPick)

	def calcVpVs(self):
		# Insert code to calc best fitting line to p travel and s-p times
		# Should be a straight line
		pass

	def _onPick(self, event):
		artist = event.artist
		label  = artist.get_label()
		self.uiReceiverDisplay.setText(label)

	def addPick(self, ptravel, stravel, dist, receiver):
		# Calculate sptime
		sptime = stravel - ptravel
		# Add a p traveltime and SP-time to the plot
		self.ptravels.append(ptravel)
		self.stravels.append(stravel)
		self.sptimes.append(sptime)
		self.dists.append(dist)
		self.receivers.append(receiver)

		self.plotWadati(receiver)
		self.plotDistance(receiver)

	def plotWadati(self, receiver):
		wadati_canvas = self.uiWadatiPlotMpl.canvas

		tolerance = 10
		for i in range(len(self.ptravels)):
			wadati_canvas.ax.scatter(self.ptravels[i], self.sptimes[i], 12, marker='o', color='k', picker=tolerance, zorder=10, label=self.receivers[i])

		wadati_canvas.draw_idle()

	def plotDistance(self, receiver):
		distance_canvas = self.uiDistancePlotMpl.canvas

		tolerance = 10
		for i in range(len(self.ptravels)):
			distance_canvas.ax.scatter(self.dists[i], self.ptravels[i], 12, marker='o', color="red", picker=tolerance, zorder=10, label=self.receivers[i])
			distance_canvas.ax.scatter(self.dists[i], self.stravels[i], 12, marker='o', color="blue", picker=tolerance, zorder=10, label=self.receivers[i])

		distance_canvas.draw_idle()

	# Function from TG - uses orthogonal distance regression
	def fitLine(self, p_times, sp_times, p_errors=None, s_errors=None):
		def func(p, x):
			x = np.array(x)
			return p[0] * x + p[1]
		linear = odr.Model(func)
		my_data = odr.RealData(p_times, sp_times, sx=p_errors, sy=sp_errors)
		my_odr  = odr.ODR(my_data, linear, beta0=[0.75, 0.])
		results = my_odr.run()


class SplittingAnalysisWindow(qt.QMainWindow):

	# ------------------------------
	# Class initialisation functions
	# ------------------------------

	def __init__(self, event, catalogue_type):
		super().__init__()

		self.event = event

		# Set up default inputs for splitting analysis
		#	np_phi 		: number of grid points for fast direction
		#	np_dt  		: number of grid points for lag time
		# 	tlag_scale  : max lag time of error surface 
		self.np_phi = 181
		if catalogue_type == "local":
			# Default maximum of 1 second delay
			self.np_dt = 101
		elif catalogue_type == "teleseismic":
			# Default maximum of 4 seconds delay
			self.np_dt = 401
		self.tlag_scale = (self.np_dt - 1) * self.event.delta

		# Create copies of the event traces
		tmp_stream = self.event.stream.copy()

		# Detrend the data
		tmp_stream = (tmp_stream.detrend("linear")).detrend("demean")

		# Grab the North and East components
		self.N = tmp_stream.select(channel="*N").data
		self.E = tmp_stream.select(channel="*E").data

		# Rotate components to Radial-Transverse
		tmp_stream.rotate(method="NE->RT", back_azimuth=self.event.baz)
		self.R = tmp_stream.select(channel="*R").data
		self.T = tmp_stream.select(channel="*T").data

		# Initialise grade variables
		self._grade = None



		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/sws_window.ui', self)

		# Populate source and receiver information boxes
		self._updateSourceInformation()
		self._updateReceiverInformation()

		self.calculateSplitting()

		### TO-DO : Plot hodograms (use sliced traces, move code into this class from sourcereceiver.py)
		###			Plot error landscape
		###			Plot corrected traces
		### 		Populate splitting results box

	# ------------------------------

	# ------------------------------
	# Connection and event functions
	# ------------------------------

	def connect(self):
		# Connect to accept/reject buttons
		self.uiAcceptButton.clicked.connect(self.accept)
		self.uiRejectButton.clicked.connect(self.reject)

		# Connect to null button
		self.uiNullButton.clicked.connect(self.setNull)

		# Connect to grading buttons
		self.uiGrade1Button.clicked.connect(lambda: self.gradeSplit("1"))
		self.uiGrade2Button.clicked.connect(lambda: self.gradeSplit("2"))
		self.uiGrade3Button.clicked.connect(lambda: self.gradeSplit("3"))
		self.uiGrade4Button.clicked.connect(lambda: self.gradeSplit("4"))

	# ------------------------------

	def rotateComponents(self):
		pass

	def calculateSplitting(self):
		### TO-DO : Grab event traces
		###		  	Rotate onto Radial and Transverse
		###			Detrend data
		###			Based on S pick, find dominant freq
		###			Create window and slice traces
		###			Run splitting analysis (c-compiled python package)
		###				Run over window start
		###					Run over window end
		###						Run over ranges of phi and dt
		###				Return all measurements of phi and dt for window sizes
		###				Return error landscape
		###			# Calculate source polarisation
		pass

	def accept(self):
		pass

	def reject(self):
		pass

	def setNull(self):
		pass

	def gradeSplit(self, grade):
		self.grade = grade
		pass

	# ------------------
	# Plotting functions
	# ------------------

	def plotTraces(self):
		# Grab a copy of the traces and rotate onto R-T
		tmp_stream = self.event.stream.copy()
		tmp_stream.rotate(method="NE->RT", back_azimuth=self.event.baz)

		R = tmp_stream.select(channel="*R")
		T = tmp_stream.select(channel="*T")

		self.uiRadialMpl.plot(R.times(type="relative"), R.data, linewidth=0.5)
		self.uiTransverseMpl.plot(T.times(type="relative"), T.data, linewidth=0.5)

		### TO-DO: Add corrected traces

	def plotClusterAnalysis(self):
		pass

	def plotHodograms(self):
		# Plot particle motion of uncorrected waveforms
		self.event.plotHodogram(self.uiHodogram1Mpl.canvas.ax)

		# Plot particle motion of corrected waveforms
		self.event.plotHodogram(self.uiHodogram2Mpl.canvas.ax, phi=self.fast_direction, dt=self.delay_time)

	def plotConfidencePlot(self):
		pass

	# ------------------

	# ---------------------------
	# Information boxes functions
	# ---------------------------

	def _updateSourceInformation(self):
		self.uiReceiverNameDisplay.setText(self.event.receiver.station)
		self.uiReceiverLonDisplay.setText(f"{self.event.receiver.longitude:.4f}")
		self.uiReceiverLatDisplay.setText(f"{self.event.receiver.latitude:.4f}")
		self.uiReceiverElevDisplay.setText(f"{self.event.receiver.elevation:.4f}")
		self.uiReceiverDepDisplay.setText(self.event.receiver.deployment.isoformat().split("T")[0])
		self.uiReceiverRetDisplay.setText(self.event.receiver.retrieval.isoformat().split("T")[0])

	def _updateReceiverInformation(self):
		self.uiOriginDateDisplay.setText(self.event.source.otime.isoformat().split("T")[0])
		self.uiOriginTimeDisplay.setText(self.event.source.otime.isoformat().split("T")[1])
		self.uiSourceLonDisplay.setText(f"{self.event.source.longitude:.4f}")
		self.uiSourceLatDisplay.setText(f"{self.event.source.latitude:.4f}")
		self.uiSourceDepthDisplay.setText(f"{self.source.depth:.4f}")
		if (type(self.event.source.magnitude) == float) or (type(self.event.source.magnitude) == np.float64):
			self.uiSourceMagDisplay.setText(f"{self.event.source.magnitude:.2f}")
		else:
			self.uiSourceMagDisplay.setText(str(self.event.source.magnitude))
		self.uiSourceIDDisplay.setText(str(self.event.source.sourceid))

	def _updateSplittingResults(self):
		self.uiFastDirDisplay.setText(self.fast_direction)
		self.uiFastDirErrDisplay.setText(self.fast_direction_error)
		self.uiDelayTimeDisplay.setText(self.delay_time)
		self.uiDelayTimeErrDisplay.setText(self.delay_time_error)
		self.uiSourcePolDisplay.setText(self.source_polarisation)

	@property
	def grade(self):
		return self._grade

	@grade.setter
	def grade(self, value):
		self._grade = value
	
	# ----------


class PickingWindow(qt.QMainWindow):

	# ------------------------------
	# Class initialisation functions
	# ------------------------------

	def __init__(self, catalogue, source=None, receiver=None, filt=None):
		super().__init__()

		self.catalogue = catalogue
		self.catalogue_path = pathlib.Path(self.catalogue.catalogue_path)
		self.source   = source
		self.receiver = receiver
		self.filt = filt
		self.default_filt = filt

		# Initialise attribute defaults
		self.pick_type = "P"
		self.pick_line_color = "red"
		self.lines = {}
		self.components = "ZNE"
		self.show_spectrogram = False

		self.sr_pairs = []
		self.counter = 0

		# If picking all receivers for a given source
		if receiver == None and source != None:
			rec_path = self.catalogue_path / "data"
			receivers = list(rec_path.glob("*/source.{}.*.z".format(source.sourceid)))
			for receiver in receivers:
				rec = receiver.parts[-1].split(".")[2]
				rec = self.catalogue.network.lookupReceiver(rec)
				rec = psm.Receiver(rec)
				self.sr_pairs.append([source, rec])

			self.wadatiWindow = WadatiWindow(self)

		# If picking all sources at a given receiver
		elif receiver != None and source == None:
			src_path = self.catalogue_path / "data" / receiver.station
			sources = list(src_path.glob("*.z"))
			for source in sources:
				src = source.parts[-1].split(".")[1]
				src = self.catalogue.lookupSource(src)
				src = psm.Source(src)
				self.sr_pairs.append([src, receiver])

		elif receiver != None and source != None:
			self.sr_pairs.append([source, receiver])

		sr = self.sr_pairs[0]
		self.current_sr = psm.SourceReceiverPair(self.catalogue_path, sr[0], sr[1])

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/trace_window.ui', self)
		self.uiToggleSpectrogramStacked.setCurrentIndex(0)

		self._updateSourceInformation(self.current_sr.source)
		self._updateReceiverInformation(self.current_sr.receiver)

		self.uiComponent1Mpl.canvas.tracePlot()
		self.uiComponent2Mpl.canvas.tracePlot()
		self.uiComponent3Mpl.canvas.tracePlot()
		self.uiSpectrogramPlot.canvas.tracePlot()

		if self.default_filt != None:
			self._populateFilter()

		self.plotTraces()

		self.connect()

		self.setWindowTitle('SeisPick - Trace picking window')
		self.setWindowIcon(QtGui.QIcon("misc/icon.png"))
		self.show()
	
	def connect(self):
		"""
		Connect to user interface callbacks.

		"""

		self.uiSaveAction.triggered.connect(self.saveTrace)

		self.uiDefaultFilterAction.triggered.connect(self.defaultFilter)
		self.uiPPickAction.triggered.connect(lambda: self.updatePick("P"))
		self.uiSPickAction.triggered.connect(lambda: self.updatePick("S"))
		self.uiCustomPickAction.triggered.connect(lambda: self.updatePick("C"))

		self.uiApplyFilterButton.clicked.connect(self.applyFilter)
		self.uiRemoveFilterButton.clicked.connect(self.removeFilter)

		self.uiRotateCompButton.clicked.connect(self.rotateComponents)
		self.uiViewSpectrogramButton.clicked.connect(self.toggleSpectrogram)
		self.uiNextTraceButton.clicked.connect(self.nextTrace)
		self.uiLastTraceButton.clicked.connect(self.previousTrace)
		self.uiRejectTraceButton.clicked.connect(self.rejectTrace)
		self.uiToggleLimitsButton.clicked.connect(self.toggleLims)
		self.uiResetPlotButton.clicked.connect(self.resetPlot)
	
	def plotconnect(self):
		"""
		Connect to plot callbacks.

		"""

		self.c1_moveid = self.uiComponent1Mpl.canvas.mpl_connect('motion_notify_event', self._onMove)
		self.c2_moveid = self.uiComponent2Mpl.canvas.mpl_connect('motion_notify_event', self._onMove)
		self.c3_moveid = self.uiComponent3Mpl.canvas.mpl_connect('motion_notify_event', self._onMove)

		self.c1_clickid = self.uiComponent1Mpl.canvas.mpl_connect('button_press_event', self._onClick)
		self.c2_clickid = self.uiComponent2Mpl.canvas.mpl_connect('button_press_event', self._onClick)
		self.c3_clickid = self.uiComponent3Mpl.canvas.mpl_connect('button_press_event', self._onClick)

		self.c1_releaseid = self.uiComponent1Mpl.canvas.mpl_connect('button_release_event', self._onRelease)
		self.c2_releaseid = self.uiComponent2Mpl.canvas.mpl_connect('button_release_event', self._onRelease)
		self.c3_releaseid = self.uiComponent3Mpl.canvas.mpl_connect('button_release_event', self._onRelease)

	def plotdisconnect(self):
		"""
		Temporarily disconnect current plot callbacks to avoid recursion.

		"""

		self.uiComponent1Mpl.canvas.mpl_disconnect(self.c1_moveid)
		self.uiComponent2Mpl.canvas.mpl_disconnect(self.c2_moveid)
		self.uiComponent3Mpl.canvas.mpl_disconnect(self.c3_moveid)

		self.uiComponent1Mpl.canvas.mpl_disconnect(self.c1_clickid)
		self.uiComponent2Mpl.canvas.mpl_disconnect(self.c2_clickid)
		self.uiComponent3Mpl.canvas.mpl_disconnect(self.c3_clickid)

		self.uiComponent1Mpl.canvas.mpl_disconnect(self.c1_releaseid)
		self.uiComponent2Mpl.canvas.mpl_disconnect(self.c2_releaseid)
		self.uiComponent3Mpl.canvas.mpl_disconnect(self.c3_releaseid)
	
	def saveTrace(self):
		"""
		Save current pick information to current source-receiver pair.

		"""

		if self.current_sr != None:
			print("Saving source...")
			self.current_sr.saveData()
		else:
			return

	def keyPressEvent(self, event):
		"""
		Override default callback for KeyPressEvent

		Parameters
		----------
		event : KeyPressEvent
			Contains information about the key press event

		"""

		if event.key() == Qt.Key_U or event.key() == Qt.Key_D:
			if not self.pick_line:
				return
			else:
				pick_polarity = chr(event.key())
				self.uiPolarityDisplay.setText(pick_polarity)
				self.current_sr.addPick(info_type="polarity", 
								 value=pick_polarity, 
								 pick_type=self.pick_type)

	def _onMove(self, event):
		"""
		Callback for motion notify event

		Parameters
		----------
		event : motion event
			Contains information about the mouse motion event

		"""

		c1_canvas = self.uiComponent1Mpl.canvas
		c2_canvas = self.uiComponent2Mpl.canvas
		c3_canvas = self.uiComponent3Mpl.canvas

		if not event.inaxes in [c1_canvas.ax, c2_canvas.ax, c3_canvas.ax]:
			return

		c1_canvas.restore_region(self.c1_bg)
		c2_canvas.restore_region(self.c2_bg)
		c3_canvas.restore_region(self.c3_bg)

		if event.key == "shift":
			if not self.trace_drag:
				return

			xpress, ypress = self.trace_click
			xmove, ymove = event.xdata, event.ydata
			dx = xmove - xpress
			dy = ymove - ypress

			rect1 = Rectangle((xpress, ypress), dx, dy, edgecolor='red', fill=False)
			rect2 = Rectangle((xpress, ypress), dx, dy, edgecolor='red', fill=False)
			rect3 = Rectangle((xpress, ypress), dx, dy, edgecolor='red', fill=False)

			c1_canvas.ax.add_patch(rect1)
			c2_canvas.ax.add_patch(rect2)
			c3_canvas.ax.add_patch(rect3)

			c1_canvas.ax.draw_artist(rect1)
			c2_canvas.ax.draw_artist(rect2)
			c3_canvas.ax.draw_artist(rect3)

		else:
			x = event.xdata

			self.c1_cursor.set_data([x, x], self.c1_y)
			self.c2_cursor.set_data([x, x], self.c2_y)
			self.c3_cursor.set_data([x, x], self.c3_y)

			c1_canvas.ax.draw_artist(self.c1_cursor)
			c2_canvas.ax.draw_artist(self.c2_cursor)
			c3_canvas.ax.draw_artist(self.c3_cursor)

		self._replotLines()				

		c1_canvas.blit(c1_canvas.ax.bbox)
		c2_canvas.blit(c2_canvas.ax.bbox)
		c3_canvas.blit(c3_canvas.ax.bbox)

	def _onClick(self, event):
		"""
		Callback for mouse click event

		Parameters
		----------
		event : click event
			Contains information about the mouse click event

		"""

		if event.key == "shift":
			if self.trace_drag:
				return
			self.trace_drag = True
			self.trace_click = event.xdata, event.ydata
			return

		c1_canvas = self.uiComponent1Mpl.canvas
		c2_canvas = self.uiComponent2Mpl.canvas
		c3_canvas = self.uiComponent3Mpl.canvas

		adjusted_xdata = round(event.xdata / self.current_sr.delta) * self.current_sr.delta

		if event.key == "control":
			if event.button == 1:
				color = "green"
				pick_type = "wb"
				self.current_sr.addPick(info_type="window", value=adjusted_xdata, pick_type="wb")
			elif event.button == 3:
				color = "green"
				pick_type = "we"
				self.current_sr.addPick(info_type="window", value=adjusted_xdata, pick_type="we")
		else:
			if event.button == 1:
				color = self.pick_line_color
				pick_type = self.pick_type
				self.current_sr.addPick(info_type="pick", value=adjusted_xdata, pick_type=self.pick_type)

				pick_time = self.current_sr.starttime + adjusted_xdata
				self.uiPickTimeDisplay.setText(pick_time.isoformat())
				self.uiPhaseDisplay.setText(self.pick_type)
			elif event.button == 3:
				# Implement pick delete option
				pass

		c1_line = c1_canvas.ax.axvline(adjusted_xdata, linewidth=1, color=color, animated=True)
		c2_line = c2_canvas.ax.axvline(adjusted_xdata, linewidth=1, color=color, animated=True)
		c3_line = c3_canvas.ax.axvline(adjusted_xdata, linewidth=1, color=color, animated=True)

		self.lines[pick_type] = [c1_line, c2_line, c3_line]

		c1_canvas.restore_region(self.c1_bg)
		c2_canvas.restore_region(self.c2_bg)
		c3_canvas.restore_region(self.c3_bg)

		self._replotLines()

		c1_canvas.blit(c1_canvas.ax.bbox)
		c2_canvas.blit(c2_canvas.ax.bbox)
		c3_canvas.blit(c3_canvas.ax.bbox)

	def _onRelease(self, event):
		c1_canvas = self.uiComponent1Mpl.canvas
		c2_canvas = self.uiComponent2Mpl.canvas
		c3_canvas = self.uiComponent3Mpl.canvas

		if not event.inaxes in [c1_canvas.ax, c2_canvas.ax, c3_canvas.ax]:
			return
		if event.key != "shift":
			return
		if not self.trace_drag:
			return

		xpress, ypress = self.trace_click
		xrelease, yrelease = event.xdata, event.ydata

		xlims = (min(xpress, xrelease), max(xpress, xrelease))
		ylims = (min(ypress, yrelease), max(ypress, yrelease))

		self.previous_lims = self.lims
		self.lims = [xlims, ylims]

		self.plotTraces(replot=True)

	def plotTraces(self, replot=False):
		c1_canvas = self.uiComponent1Mpl.canvas
		c2_canvas = self.uiComponent2Mpl.canvas
		c3_canvas = self.uiComponent3Mpl.canvas

		c1_canvas.ax.clear()
		c2_canvas.ax.clear()
		c3_canvas.ax.clear()

		# Create background variables
		self.c1_bg = None
		self.c2_bg = None
		self.c3_bg = None

		self.trace_drag  = False
		self.trace_click = None

		if not replot:
			# Add limits variable
			self.previous_lims = None
			self.lims = None

			# Set filter to default filter
			self.filt = self.default_filt

			# Create an instance of the Source class
			self.current_sr = psm.SourceReceiverPair(self.catalogue_path, self.source, self.receiver)

			# Look up any picks and make vertical lines
			for pick, data in self.current_sr.picks.items():
				pick_time = data["rtime"]
				if pick == "we" or pick == "wb":
					color = "green"
				else:
					color = "gray"

				c1_canvas.ax.axvline(pick_time, linewidth=0.75, color=color)
				c2_canvas.ax.axvline(pick_time, linewidth=0.75, color=color)
				c3_canvas.ax.axvline(pick_time, linewidth=0.75, color=color)

		else:
			self._replotLines()

		self.current_sr.filter(filt=self.filt)

		# Plot the traces
		self.current_sr.plotTraces([c1_canvas.ax, c2_canvas.ax, c3_canvas.ax], lims=self.lims)

		# Connect to each trace to grab the background once Qt has done resizing
		c1_canvas.mpl_connect('draw_event', self._c1DrawEvent)
		c2_canvas.mpl_connect('draw_event', self._c2DrawEvent)
		c3_canvas.mpl_connect('draw_event', self._c3DrawEvent)

		# Draw when Qt has done all resizing
		c1_canvas.draw_idle()
		c2_canvas.draw_idle()
		c3_canvas.draw_idle()

		# Initialise the cursor to track mouse position on the axes
		self.c1_cursor = c1_canvas.ax.axvline(5, linewidth=1, color='0.5', animated=True)
		_, self.c1_y = self.c1_cursor.get_data()
		self.c2_cursor = c2_canvas.ax.axvline(5, linewidth=1, color='0.5', animated=True)
		_, self.c2_y = self.c2_cursor.get_data()
		self.c3_cursor = c3_canvas.ax.axvline(5, linewidth=1, color='0.5', animated=True)
		_, self.c3_y = self.c3_cursor.get_data()

		self.plotconnect()

	def _c1DrawEvent(self, event):
		"""
		Grab the trace background when it is drawn

		Parameters
		----------
		event : DrawEvent
			Contains information about the draw event

		"""

		self.c1_bg = event.canvas.copy_from_bbox(event.canvas.ax.bbox)

	def _c2DrawEvent(self, event):
		"""
		Grab the trace background when it is drawn

		Parameters
		----------
		event : DrawEvent
			Contains information about the draw event

		"""

		self.c2_bg = event.canvas.copy_from_bbox(event.canvas.ax.bbox)

	def _c3DrawEvent(self, event):
		"""
		Grab the trace background when it is drawn

		Parameters
		----------
		event : DrawEvent
			Contains information about the draw event

		"""

		self.c3_bg = event.canvas.copy_from_bbox(event.canvas.ax.bbox)

	def _replotLines(self):
		c1_canvas = self.uiComponent1Mpl.canvas
		c2_canvas = self.uiComponent2Mpl.canvas
		c3_canvas = self.uiComponent3Mpl.canvas

		for pick_type, lines in self.lines.items():
			c1_canvas.ax.draw_artist(lines[0])
			c2_canvas.ax.draw_artist(lines[1])
			c3_canvas.ax.draw_artist(lines[2])

	def _updateReceiverInformation(self, receiver):
		"""
		Populates the Receiver Info widget with receiver information

		"""

		self.uiReceiverNameDisplay.setText(receiver.station)
		self.uiReceiverLonDisplay.setText(f"{receiver.longitude:.4f}")
		self.uiReceiverLatDisplay.setText(f"{receiver.latitude:.4f}")
		self.uiReceiverElevDisplay.setText(f"{receiver.elevation:.4f}")
		self.uiReceiverDepDisplay.setText(receiver.deployment.isoformat().split("T")[0])
		self.uiReceiverRetDisplay.setText(receiver.retrieval.isoformat().split("T")[0])

	def _updateSourceInformation(self, source):
		"""
		Populates the Source Info widget with source information

		"""

		self.uiOriginDateDisplay.setText(source.otime.isoformat().split("T")[0])
		self.uiOriginTimeDisplay.setText(source.otime.isoformat().split("T")[1])
		self.uiSourceLonDisplay.setText(f"{source.longitude:.4f}")
		self.uiSourceLatDisplay.setText(f"{source.latitude:.4f}")
		self.uiSourceDepthDisplay.setText(f"{source.depth:.4f}")
		if (type(source.magnitude) == float) or (type(source.magnitude) == np.float64):
			self.uiSourceMagDisplay.setText(f"{source.magnitude:.2f}")
		else:
			self.uiSourceMagDisplay.setText(str(source.magnitude))
		self.uiSourceIDDisplay.setText(str(source.sourceid))

	def defaultFilter(self):
		"""
		Opens a default filter window.

		"""

		self.defaultFilterDialogue = DefaultFilterDialogue(self)

		if not self.defaultFilterDialogue.exec_():
			return

	def applyFilter(self):
		"""
		Applies the 

		TO-DO
		-----
		Parse for each filter type (bandpass, lowpass, highpass etc)

		"""

		try:
			filt_type = self.uiFilterTypeComboBox.currentText()
			corners  = int(self.uiPoleCountComboBox.currentText())
			zerophase = self.uiZeroPhaseCheckBox.isChecked()
			freqmin  = float(self.uiLowFreqInput.text())
			freqmax = float(self.uiHighFreqInput.text())

			self.filt = {'type': filt_type,
						 'corners': corners,
						 'zerophase': zerophase,
						 'freqmin': freqmin,
						 'freqmax': freqmax}
		except ValueError:
			qt.QMessageBox.about(self, "Warning!", "You appear to have tried to use an incomplete/incorrect filter. Fill in all of the options and try again.")
			return

		self.plotdisconnect()
		self.plotTraces(replot=True)

	def removeFilter(self):
		"""
		Remove any filters applied to the current source-receiver pair.

		"""

		self.plotdisconnect()
		self.filt = None
		self.current_sr.filter()
		self.plotTraces(replot=True)

	def _populateFilter(self):
		"""
		Populate filter input with filter parameters.

		"""

		filter_type = self.default_filt["type"]

		self.uiFilterTypeComboBox.setCurrentText(filter_type)
		self.uiPoleCountComboBox.setCurrentText(str(self.default_filt["corners"]))

		if filter_options["zerophase"]:
			self.uiZeroPhaseCheckBox.setChecked(True)

		if filter_type == "bandpass":
			self.uiLowFreqInput.setText(str(self.default_filt["minfreq"]))
			self.uiHighFreqInput.setText(str(self.default_filt["maxfreq"]))

		if filter_type == "lowpass":
			self.uiLowFreqInput.setText(str(self.default_filt["freq"]))

		if filter_type == "highpass":
			self.uiHighFreqInput.setText(str(self.default_filt["freq"]))

	def updatePick(self, pick_type):
		"""
		Updates pick type and corresponding variables

		Parameters
		----------
		pick_type : str
			String specifying the type of pick ("S", "P", "C")

		"""

		if pick_type == "P":
			self.pick_type = pick_type
			self.pick_line_color = "red"

		elif pick_type == "S":
			self.pick_type = pick_type
			self.pick_line_color = "blue"

		elif pick_type == "C":
			# Only run if the radio button is being toggled on
			if not self.uiCustomRadio.isChecked():
				return
			self.customPickDialogue = CustomPickDialogue()

			if not self.customPickDialogue.exec_():
				return

		self.uiPolarityDisplay.setText("")

	def rotateComponents(self):
		"""
		Rotate source-receiver components.

		"""

		if self.components == "ZNE":
			self.components = "ZRT"
			self.current_sr.rotate(method="NE->RT")

		elif self.components == "ZRT":
			self.components = "ZNE"
			self.current_sr.rotate(method="RT->NE")

		self.plotTraces(replot=True)

	def toggleSpectrogram(self):
		"""
		Toggles spectrogram view

		"""

		self.show_spectrogram = not self.show_spectrogram
		self.plotSpectrogram()

	def plotSpectrogram(self):
		"""
		Plot spectrogram for current event

		"""

		if self.show_spectrogram:
			self.uiToggleSpectrogramStacked.setCurrentIndex(1)
			self.current_sr.plotSpectrogram(self.uiSpectrogramPlot.canvas.ax)

		else:
			self.uiToggleSpectrogramStacked.setCurrentIndex(0)

	def nextTrace(self, trace_removed=False):
		"""
		Moves the counter forward and updates the current source-receiver pair.

		"""

		if not trace_removed:
			self.counter += 1
		else:
			# Add test to see if source-receiver list is empty.
			pass

		if self.counter == len(self.sr_pairs):
			self.counter = 0

		sr = self.sr_pairs[self.counter]
		self.current_sr = psm.SourceReceiverPair(self.catalogue_path, sr[0], sr[1])

		self._updateSourceInformation(sr[0])
		self._updateReceiverInformation(sr[1])

		self.updateTrace()

	def previousTrace(self):
		"""
		Moves the counter back and updates the current source-receiver pair.

		"""

		self.counter += -1

		if self.counter == -1:
			self.counter = len(self.sr_pairs) - 1

		sr = self.sr_pairs[self.counter]
		self.current_sr = psm.SourceReceiverPair(self.catalogue_path, sr[0], sr[1])

		self._updateSourceInformation(sr[0])
		self._updateReceiverInformation(sr[1])

		self.updateTrace()

	def updateTrace(self):
		# Check if P and S have been picked (and window exists)
		if "P_manual" in self.current_sr.picks.keys() and "S_manual" in self.current_sr.picks.keys() and self.wadatiWindow:
			ptravel = self.current_sr.starttime - self.current_sr.otime
			ptravel += self.current_sr.picks["P_manual"]["rtime"]
			stravel = self.current_sr.starttime - self.current_sr.otime
			stravel += self.current_sr.picks["S_manual"]["rtime"]
			self.wadatiWindow.addPick(ptravel, stravel, self.receiver)

		self.saveTrace()

		self.uiPRadio.setChecked(True)

		self.lines = {}
		self.pick_time = ""
		self.pick_phase = ""

		self.uiPolarityDisplay.setText("")
		self.uiPickTimeDisplay.setText(self.pick_time)
		self.uiPhaseDisplay.setText(self.pick_phase)

		self.plotdisconnect()
		self.plotTraces()

	def rejectTrace(self):
		"""
		Reject the current trace.

		Unlinks the source-receiver data files and removes the arrival
		from the catalogue.

		"""

		if self.current_sr != None:
			print("Deleting trace...")
			self.current_sr.removeFiles()
			self.catalogue.removeArrival(self.current_sr.source.sourceid, self.current_sr.receiver.receiverid)

			del self.sr_pairs[self.counter]

			self.nextTrace(trace_removed=True)

	def toggleLims(self):
		"""
		Toggles plot limits

		TO-DO
		-----
		Implement a proper system of undo/redo

		"""

		self.lims, self.previous_lims = self.previous_lims, self.lims

		self.plotTraces(replot=True)

	def resetPlot(self):
		"""
		Reset the plot of the traces, returning to the default limits and
		filter.

		"""

		self.plotdisconnect()
		self.lims = None
		self.filt = self.default_filter
		self.plotTraces()


if __name__ == "__main__":

	app = qt.QApplication(sys.argv)
	seispick = SeisPick()
	sys.exit(app.exec_())