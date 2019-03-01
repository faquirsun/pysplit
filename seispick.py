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

import metainfo as psm
import catalogue as cat

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

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

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/main.ui', self)

		# Initialise toggle trackers
		self.ctrl  = False
		self.shift = False

		# Initialise empty dictionary for catalogue parameters
		self.catalogue_parameters = {}

		# Connect to actions and interactive widgets
		self.connect()

		# Hide the map widget initially
		self.uiMapMpl.hide()
		self.uiStatusStacked.setCurrentIndex(1)

		# Set initial filter to None
		self.filt = None

		self.setWindowTitle('SeisPick - Main catalogue window')
		self.setWindowIcon(QtGui.QIcon("misc/icon.png"))
		self.show()

	# ------------------------------

	# ------------------------------
	# Connection and event functions
	# ------------------------------

	def connect(self):
		# File menu connections
		self.uiNewCatalogueAction.triggered.connect(self.newCatalogue)
		self.uiLoadCatalogueAction.triggered.connect(self.loadCatalogue)

		# Settings menu connections
		self.uiDefaultFilterAction.triggered.connect(self.defaultFilter)

		# Map connections
		self.uiMapMpl.canvas.mpl_connect('pick_event', self._onPick)
		self.uiMapMpl.canvas.mpl_connect('button_press_event', self._onMapClick)
		self.uiMapMpl.canvas.mpl_connect('motion_notify_event', self._onMapMove)
		self.uiMapMpl.canvas.mpl_connect('button_release_event', self._onMapRelease)

		# Receiver list connection
		self.uiReceiversListView.doubleClicked.connect(self.receiverSelect)

		# Source list connection
		self.uiSourcesListView.doubleClicked.connect(self.sourceSelect)

		# Catalogue load button connection
		self.uiReceiverLoadButton.clicked.connect(self.loadReceivers)
		self.uiSourceLoadButton.clicked.connect(self.loadSources)
		self.uiArrivalLoadButton.clicked.connect(self.loadArrivals)
		self.uiWaveformLoadButton.clicked.connect(self.loadWaveforms)

		# Plot option tickboxes connection
		self.uiPickedCheckBox.stateChanged.connect(self.pickStateChanged)
		self.uiUnpickedCheckBox.stateChanged.connect(self.unpickedStateChanged)
		self.uiReceiversCheckBox.stateChanged.connect(self.receiverStateChanged)

		self.uiPlotMapButton.clicked.connect(lambda: self.plotCatalogueMap(replot=True))
		self.uiResetMapButton.clicked.connect(self.plotCatalogueMap)

	def keyPressEvent(self, event):
		# Toggle Ctrl Modifier on
		if event.key() == Qt.Key_Control:
			self.ctrl = True

		# Toggle Shift Modifier on
		if event.key() == Qt.Key_Shift:
			self.shift = True

	def keyReleaseEvent(self, event):
		# Toggle Ctrl Modifier off
		if event.key() == Qt.Key_Control:
			self.ctrl = False

		# Toggle Shift Modifier off
		if event.key() == Qt.Key_Shift:
			self.shift = False

	def defaultFilter(self):
		self.defaultFilterDialogue = DefaultFilterDialogue(self)

		if not self.defaultFilterDialogue.exec_():
			return

	# ------------------------------

	# ---------------------------------
	# Catalogue instantiation functions
	# ---------------------------------

	def newCatalogue(self):
		self.newCatalogueDialogue = NewCatalogueDialogue(self)

		if not self.newCatalogueDialogue.exec_():
			return

		self.createCatalogue(new=True)

	def loadCatalogue(self):
		# Open file dialogue and choose catalogue directory
		self.catalogue_path = qt.QFileDialog.getExistingDirectory(self, 'Choose catalogue directory')

		# Test validity of the catalogue path
		if not (os.path.exists("{}/data".format(self.catalogue_path))) and (os.path.exists("{}/metafiles".format(self.catalogue_path))):
			qt.QMessageBox.about(self, "Error!", "The directory you have chosen is invalid, please try again.")
			return

		# Parse the catalogue metafile
		try:
			self._parseCatalogueMetafile()
		except FileNotFoundError:
			qt.QMessageBox.about(self, "Error!", "Unable to find file containing the metadata for this catalogue. Please double-check it exists.")
			return

		self.createCatalogue()

	def createCatalogue(self, new=False):
		# Handles catalogue creation and interface population
		# XXX
		# Add function to reset all variables so as to avoid clashes
		# Also clear all Display objects etc
		# XXX

		# Display catalogue information
		self.uiCatalogueTypeDisplay.setText(self.catalogue_parameters["catalogue_type"])
		self.uiCatalogueNameDisplay.setText(self.catalogue_parameters["catalogue_name"])
		self.uiStatusStacked.setCurrentIndex(0)

		# Enable the receiver and source selection lists
		self.uiReceiversListView.setEnabled(True)
		self.uiSourcesListView.setEnabled(True)

		if self.catalogue_parameters["catalogue_type"] == "teleseismic":
			self.uiStatusBar.showMessage("Creating an instance of TeleseismicCatalogue...")
			# Display catalogue information
			self.uiPlotOptionsStacked.setCurrentIndex(1)
			# Create an instance of the Teleseismic Catalogue class
			self.catalogue = cat.TeleseismicCatalogue(new, **self.catalogue_parameters)

		elif self.catalogue_parameters["catalogue_type"] == "local":
			self.uiStatusBar.showMessage("Creating an instance of LocalCatalogue...")
			# Display local catalogue progress page
			self.uiPlotOptionsStacked.setCurrentIndex(0)
			# Create an instance of the Local Catalogue class
			self.catalogue = cat.LocalCatalogue(new, **self.catalogue_parameters)

		if new:
			# Populate the new catalogue with parameters
			self.catalogue.newCatalogue(self.catalogue_parameters)
			self.uiStatusBar.showMessage("Generating new catalogue...")
			# Generate the catalogue
			self.catalogue.generateCatalogue()
		else:
			pass
			# Load catalogue with parameters
			# self.catalogue.loadCatalogue(self.catalogue_parameters)
			# self.uiStatusBar.showMessage("Loading existing catalogue files...")
			# self.loadReceivers()
			# self.loadSources()

		# Populate the information wiedgets
		self._populateReceiverList()
		if self.catalogue_parameters["catalogue_type"] == "teleseismic":
			self._populateCatalogueInformation()

		# Plot the map of the catalogue
		self.plotCatalogueMap()

	def loadSources(self):
		sources = self.catalogue.loadSources()
		if sources:
			self.uiSourceLoadButton.setEnabled(False)
		else:
			qt.QMessageBox.about(self, "Error!", "Error loading source file.")
			return

	def loadReceivers(self):
		receivers = self.catalogue.loadReceivers()
		if receivers:
			self.uiReceiverLoadButton.setEnabled(False)
		else:
			qt.QMessageBox.about(self, "Error!", "Error loading receiver file.")
			return

	def loadArrivals(self):
		arrivals = self.catalogue.loadArrivals()
		if arrivals:
			self.uiArrivalLoadButton.setEnabled(False)
			self.uiWaveformLoadButton.setEnabled(True)
		else:
			self.uiStatusBar.showMessage("The arrivals file does not exist - attempting to generate one...")
			if self.catalogue.catalogue_type == "local":
				# Load the arrivals from the local input file provided
				self.catalogue.getArrivals(input_file=self.local_input, input_type=self.data_source)

			elif self.catalogue.catalogue_type == "teleseismic":
				self.input_phases = []
				
				# Open pop-up to get the phases to be collected
				self.telePhaseDialogue = TelePhaseDialogue(self)

				if not self.telePhaseDialogue.exec_():
					return

				self.catalogue.getArrivals(phases=self.input_phases)
				self._populateReceiverList()

			self.loadArrivals()

	def loadWaveforms(self):
		waveforms = self.catalogue.loadWaveforms()
		if waveforms:
			self.uiWaveformLoadButton.setEnabled(False)
		else:
			self.catalogue.getWaveforms()
			#qt.QMessageBox.about(self, "Error!", "Some error in there somewhere yo")
			return

		self.uiWaveformStatusDisplay.setText("Waveforms generated.")
		self.uiWaveformStatusDisplay.setStyleSheet("color: rgb(0, 255, 0)")
		self.uiWaveformLoadButton.setEnabled(False)

	def _populateReceiverList(self):
		# Populate the receiver list
		model = QtGui.QStandardItemModel(self.uiReceiversListView)
		for receiver in self.catalogue.network.receivers.name.values:
			item = QtGui.QStandardItem(receiver)
			model.appendRow(item)
		self.uiReceiversListView.setModel(model)

	def _populateCatalogueInformation(self):
		self.uiStartDateDisplay.setText(self.catalogue_parameters["start_date"])
		self.uiEndDateDisplay.setText(self.catalogue_parameters["end_date"])
		self.uiCenLonDisplay.setText(self.catalogue_parameters["clon"])
		self.uiCenLatDisplay.setText(self.catalogue_parameters["clat"])
		self.uiMinRadDisplay.setText(self.catalogue_parameters["minrad"])
		self.uiMaxRadDisplay.setText(self.catalogue_parameters["maxrad"])
		self.uiMinMagDisplay.setText(self.catalogue_parameters["minmag"])

		# Calculate number of sources and update label
		no_sources = str(len(self.catalogue.src_df.index))
		self.uiSourceCountDisplay.setText(no_sources)

	def _parseCatalogueMetafile(self):
		# Update status bar
		self.uiStatusBar.showMessage('Parsing catalogue metafile...')

		# Empty list for collecting parameters
		params = []

		# Read in the parameters from the file
		for line in open("{}/metafiles/catalogue_metafile.txt".format(self.catalogue_path), 'r'):
			if "?" in line:
				line = line.split(" ? ")
				if line[1].endswith("\n"):
					param = line[1][:-1]
				params.append(param)


		catalogue_parameters = {'catalogue_name': params[0],
								'catalogue_type': params[1],
								'catalogue_path': params[2],
								'data_source': params[3],
								'cdate': params[4],
								'archive_path': params[5],
								'archive_type': params[6],
								'archive_format': params[7],
								'rec_file': params[8],
								'start_date': params[9],
								'end_date':params[10]}

		if catalogue_parameters['catalogue_type'] == "local":
			d_specific = {'local_input': params[11]}

		if catalogue_parameters['catalogue_type'] == "teleseismic":
			d_specific = {'minmag': params[11],
			 			  'clon': params[12],
					 	  'clat': params[13],
						  'minrad': params[14],
						  'maxrad': params[15]}

		self.catalogue_parameters = {**catalogue_parameters, **d_specific}

	# ---------------------------------

	# -------------------------
	# Handlers for map plotting
	# -------------------------

	def plotCatalogueMap(self, replot=False):
		# Transmit a message to the status bar
		self.uiStatusBar.showMessage("Plotting catalogue map...")

		# Create background, lock and click variables
		self.map_background = None
		self._map_drag_lock = None
		self.map_click      = None

		# Show the matplotlib widget and clear it
		self.uiMapMpl.show()

		# If replotting, grab the new parameters
		if replot:
			self.uiMapMpl.canvas.ax.clear()
			self.catalogue.plotGeographic(self.uiMapMpl,
										   lon0=float(self.uiMinLonInput.text()),
										   lon1=float(self.uiMaxLonInput.text()),
										   lat0=float(self.uiMinLatInput.text()),
										   lat1=float(self.uiMaxLatInput.text()))

			if self.uiReceiversCheckBox.isChecked():
				self.catalogue.plotReceivers(self.uiMapMpl)

		else:
			# Plot the sources
			self.catalogue.plotGeographic(self.uiMapMpl)

			# Plot the receivers
			if self.uiReceiversCheckBox.isChecked():
				self.catalogue.plotReceivers(self.uiMapMpl)

			# Grab the latitude and longitude of the map and send to input options
			if self.catalogue.catalogue_type == "local":
				self.uiMinLonInput.setText(str(f"{self.catalogue.lon0:.5f}"))
				self.uiMaxLonInput.setText(str(f"{self.catalogue.lon1:.5f}"))
				self.uiMinLatInput.setText(str(f"{self.catalogue.lat0:.5f}"))
				self.uiMaxLatInput.setText(str(f"{self.catalogue.lat1:.5f}"))

		# Connect to the map to grab background once Qt has done resizing
		self.uiMapMpl.canvas.mpl_connect('draw_event', self._mapDrawEvent)

		# Draw the canvas once Qt has done all resizing
		self.uiMapMpl.canvas.draw_idle()

		# Transmit a message to the status bar
		self.uiStatusBar.showMessage("Catalogue load complete.")

	def _mapDrawEvent(self, event):
		# Grab the map background when it is drawn
		self.map_background = event.canvas.copy_from_bbox(event.canvas.ax.bbox)

	def unpickedStateChanged(self, int):
		if self.uiUnpickedCheckBox.isChecked():
			pass
			#self.catalogue.receiver_points.set_visible(True)
		else:
			pass
			#self.catalogue.receiver_points.set_visible(False)

	def receiverStateChanged(self, int):
		# Toggle receiver display 
		self.plotCatalogueMap(replot=True)

	def pickStateChanged(self, int):
		if self.uiPickedCheckBox.isChecked():
			pass
			#self.catalogue.pick_points.set_visible(True)
		else:
			pass
			#self.catalogue.pick_points.set_visible(False)

	def _onPick(self, event):
		# Retrieve the artist and information contained therein from the mouse event
	    artist = event.artist
	    xy = artist.get_offsets()
	    label = artist.get_label()
	    ind = event.ind

	    print(event.mouseevent.key)

	    # Convert the x/y positions of the mouse event to lat/lon using the relevant transform
	    lon, lat = xy[ind[0]][0], xy[ind[0]][1]

	    # If the data point selected is a receiver
	    if "REC" in label:
	    	rec = label.split(": ")[1]
	    	self.receiver = psm.Receiver(self.catalogue.network.receivers.query('name == @rec'))
	    	self._updateReceiverInformation()

	    	# If in receiver pick mode, open a Picking Window
	    	if self.uiPickReceiversAction.isChecked():
	    		self.pickWindow = PickingWindow(self, self.catalogue, receiver=self.receiver)

	    # If the data point selected is an source
	    if "SOURCE" in label:
	    	src = label.split(": ")[1]
	    	self.source = psm.Source(self.catalogue.src_df.query('sourceid == @src'))
	    	self._updateSourceInformation()

	    	# If in source pick mode, open a Picking Window
	    	if self.uiPickSourcesAction.isChecked():
		    	self.pickWindow = PickingWindow(self, self.catalogue, source=self.source)

	def _onMapClick(self, event):
		print(event.key)
		if not self.shift:
			return
		if self.catalogue.catalogue_type == "teleseismic":
			return
		if self._map_drag_lock is not None:
			return
		if event.inaxes != self.uiMapMpl.canvas.ax:
			return

		# Grab the x and y positions of the click point
		xpress, ypress = event.xdata, event.ydata
		self._map_drag_lock = self

		# Store the click point 
		self.map_click = xpress, ypress

	def _onMapMove(self, event):
		if not self.shift:
			return
		if self.catalogue.catalogue_type == "teleseismic":
			return
		if self._map_drag_lock is not self:
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
		self.uiMapMpl.canvas.restore_region(self.map_background)

		# Set rectangle values
		map_rectangle = Rectangle((xpress, ypress), dx, dy,
								  edgecolor='red', fill=False)
		self.uiMapMpl.canvas.ax.add_patch(map_rectangle)
		self.uiMapMpl.canvas.ax.draw_artist(map_rectangle)

		# Blit the redrawn area
		self.uiMapMpl.canvas.blit(self.uiMapMpl.canvas.ax.bbox)

	def _onMapRelease(self, event):
		if not self.shift:
			return
		if self.catalogue.catalogue_type == "teleseismic":
			return
		if self._map_drag_lock is not self:
			return

		# Grab the x and y data positions of the click point
		xpress, ypress = self.map_click

		# Reset map click and lock variables
		self.map_click      = None
		self._map_drag_lock = None

		# Grab the x and y data positions of the release point
		xrelease, yrelease = event.xdata, event.ydata
		dx = xrelease - xpress
		dy = yrelease - ypress

		# Draw background from pixel buffer
		self.uiMapMpl.canvas.restore_region(self.map_background)

		# Draw the final rectangle
		map_rectangle = Rectangle((xpress, ypress), dx, dy,
								  edgecolor='red', fill=False)
		self.uiMapMpl.canvas.ax.add_patch(map_rectangle)
		self.uiMapMpl.canvas.ax.draw_artist(map_rectangle)

		# Blit the redrawn area
		self.uiMapMpl.canvas.blit(self.uiMapMpl.canvas.ax.bbox)

		# Convert the x and y positions to lon/lat
		lonpress, latpress     = xpress, ypress
		lonrelease, latrelease = xrelease, yrelease

		# Set the text values of the lon/lat input boxes
		self.uiMinLonInput.setText(str(f"{min(xpress, xrelease):.2f}"))
		self.uiMaxLonInput.setText(str(f"{max(xpress, xrelease):.2f}"))
		self.uiMinLatInput.setText(str(f"{min(ypress, yrelease):.2f}"))
		self.uiMaxLatInput.setText(str(f"{max(ypress, yrelease):.2f}"))

	# -------------------------

	# -------------------------
	# Source/receiver functions
	# -------------------------

	def receiverSelect(self, index):
		# Parse receiver selected
		receiver = self.uiReceiversListView.model().data(index)

		# Lookup receiver information and create Receiver object
		self.receiver = psm.Receiver(self.catalogue.network.receivers.query('name == @receiver'))
		self.list_receiver = self.receiver

		# Load receiver information and print to display
		self._updateReceiverInformation()

		sources = glob.glob('{}/data/{}/*.z'.format(self.catalogue.catalogue_path, receiver.upper()))
		sources.sort(key=natural_keys)

		# Populate the source list
		model = QtGui.QStandardItemModel(self.uiSourcesListView)
		for source in sources:
			head, tail = os.path.split(source)
			item = QtGui.QStandardItem("Source {}".format(tail.split(".")[1]))
			model.appendRow(item)
		self.uiSourcesListView.setModel(model)

	def sourceSelect(self, index):
		# Parse source selected
		source = self.uiSourcesListView.model().data(index).split(" ")[1]

		# Look up source information and create Source object
		self.source = psm.Source(self.catalogue.src_df.query('sourceid == @source'))

		# Load source information and print to display
		self._updateSourceInformation()

		# Open up picking window
		self.pickWindow = PickingWindow(self, self.catalogue, receiver=self.list_receiver, source=self.source)

	def _updateReceiverInformation(self):
		self.uiReceiverNameDisplay.setText(self.receiver.station)
		self.uiReceiverLonDisplay.setText(f"{self.receiver.longitude:.4f}")
		self.uiReceiverLatDisplay.setText(f"{self.receiver.latitude:.4f}")
		self.uiReceiverElevDisplay.setText(f"{self.receiver.elevation:.4f}")
		self.uiReceiverDepDisplay.setText(self.receiver.deployment.isoformat().split("T")[0])
		self.uiReceiverRetDisplay.setText(self.receiver.retrieval.isoformat().split("T")[0])

	def _updateSourceInformation(self):
		self.uiOriginDateDisplay.setText(self.source.otime.isoformat().split("T")[0])
		self.uiOriginTimeDisplay.setText(self.source.otime.isoformat().split("T")[1])
		self.uiSourceLonDisplay.setText(f"{self.source.longitude:.4f}")
		self.uiSourceLatDisplay.setText(f"{self.source.latitude:.4f}")
		self.uiSourceDepthDisplay.setText(f"{self.source.depth:.4f}")
		if (type(self.source.magnitude) == float) or (type(self.source.magnitude) == np.float64):
			self.uiSourceMagDisplay.setText(f"{self.source.magnitude:.2f}")
		else:
			self.uiSourceMagDisplay.setText(str(self.source.magnitude))
		self.uiSourceIDDisplay.setText(str(self.source.sourceid))

	# -------------------------


class NewCatalogueDialogue(qt.QDialog):

	# ------------------------------
	# Class initialisation functions
	# ------------------------------

	def __init__(self, parent):
		super().__init__()

		self.parent = parent

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
		self.catalogue_path = qt.QFileDialog.getExistingDirectory(self, 'Choose catalogue directory')
		self.uiCatalogueDirectoryInput.setText(self.catalogue_path)

	def browseArchive(self):
		self.archive_path = qt.QFileDialog.getExistingDirectory(self, 'Choose archive directory')
		self.uiArchiveDirectoryInput.setText(self.archive_path)

	def browseReceivers(self):
		filename = qt.QFileDialog.getOpenFileName(self, 'Open file')
		self.rec_file = filename[0]
		self.uiReceiverFileInput.setText(self.rec_file)

	def browseLocalFile(self):
		filename = qt.QFileDialog.getOpenFileName(self, 'Open file')
		self.input_file = filename[0]
		self.uiLocalFileInput.setText(self.input_file)

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
			archive_format = "{year}/{jday}/*_{receiver}_{comp}*"
		
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
		self.parent.catalogue_parameters = {**catalogue_parameters, **d_specific}

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

	def __init__(self, parent):
		super().__init__()

		self.parent = parent

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

		# Bandpass filter
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
			filter_options = {'freqmin': freqmin,
							  'freqmax': freqmax,
							  'corners': corners,
							  'zerophase': zerophase}

		# Lowpass filter
		elif filter_type == "lowpass":
			try:
				freq = float(self.uiLowFreqInput.text())
			except ValueError:
				qt.QMessageBox.about(self, "Error!", "Invalid lowpass frequency.")
				return

			# Create filter options dictionary
			filter_options = {'freq': freq,
							  'corners': corners,
							  'zerophase': zerophase}

		# Highpass filter
		elif filter_type == "highpass":
			try:
				freq = float(self.uiHighFreqInput.text())
			except ValueError:
				qt.QMessageBox.about(self, "Error!", "Invalid highpass frequency.")
				return

			# Create filter options dictionary
			filter_options = {'freq': freq,
							  'corners': corners,
							  'zerophase': zerophase}

		# Update the default filter in parent class
		self.parent.filt = [filter_type, filter_options]

		# Send accept signal to Dialog
		self.accept()

	def actionReject(self):
		# Send reject signal to Dialog
		self.reject()

	# ------------------------------


class CustomPickDialogue(qt.QDialog):

	# ------------------------------
	# Class initialisation functions
	# ------------------------------

	def __init__(self, parent):
		super(CustomPickDialogue, self).__init__()

		self.parent = parent

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/custom_phase_dialogue.ui', self)

		self.setWindowTitle('SeisPick - Set custom phase pick')
		self.setWindowIcon(QtGui.QIcon("misc/icon.png"))
		self.show()

	# ------------------------------

	# ------------------------------
	# Accept/Reject action overrides
	# ------------------------------

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
		# Send reject signal to Dialog
		self.reject()

	# ------------------------------


class WadatiWindow(qt.QMainWindow):

	# ------------------------------
	# Class initialisation functions
	# ------------------------------

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

		# Set plot details (axes labels etc)
		self.uiWadatiPlotMpl.canvas.ax.set_xlabel("P traveltime / s", fontsize=10)
		self.uiWadatiPlotMpl.canvas.ax.set_ylabel("S - P traveltime / s", fontsize=10)

		# Connect to buttons and plots
		self.connect()

		self.setWindowTitle('SeisPick - Wadati plot window')
		self.setWindowIcon(QtGui.QIcon("misc/icon.png"))
		self.show()

	# ------------------------------

	# ------------------------------
	# Connection and event functions
	# ------------------------------

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

	# ------------------------------

	# --------------------------
	# Handler for new pick input
	# --------------------------

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

	# --------------------------

	# ------------------
	# Plotting functions
	# ------------------

	def plotWadati(self, receiver):
		wadati_canvas = self.uiWadatiPlotMpl.canvas

		# Clear the canvas
		wadati_canvas.ax.clear()

		# Set plot details (axes labels etc)
		wadati_canvas.ax.set_xlabel("P traveltime / s", fontsize=10)
		wadati_canvas.ax.set_ylabel("S - P traveltime / s", fontsize=10)

		# Try rescaling the image now
		wadati_canvas.ax.set_aspect('auto')

		tolerance = 10
		for i in range(len(self.ptravels)):
			wadati_canvas.ax.scatter(self.ptravels[i], self.sptimes[i], 12, marker='o', color='k', picker=tolerance, zorder=10, label=self.receivers[i])

		wadati_canvas.draw_idle()

	def plotDistance(self, receiver):
		distance_canvas = self.uiDistancePlotMpl.canvas

		# Clear the canvas
		distance_canvas.ax.clear()

		# Set plot details (axes labels etc)
		distance_canvas.ax.set_xlabel("Epicentral distance / km", fontsize=10)
		distance_canvas.ax.set_ylabel("Traveltime / s", fontsize=10)

		# Rescale the image now
		distance_canvas.ax.set_aspect('auto')

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

	# ------------------


class SplittingAnalysisWindow(qt.QMainWindow):

	# ------------------------------
	# Class initialisation functions
	# ------------------------------

	def __init__(self, event, cataloge_type):
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

		# Populate source and receiver information boxes
		self._updateSourceInformation()
		self._updateReceiverInformation()

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/sws_window.ui', self)

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

	# ---------------------------

	# ----------
	# Properties
	# ----------

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

	def __init__(self, parent, catalogue, source=None, receiver=None):
		super().__init__()

		# Distribute access to the main class variables
		self.parent = parent

		# Parse args
		self.catalogue = catalogue
		self.catalogue_path = self.catalogue.catalogue_path
		self.catalogue_name = self.catalogue.catalogue_name
		self.source   = source
		self.receiver = receiver

		# Set default filter
		self.default_filt = self.parent.filt
		self.filt = self.parent.filt

		# Set default pick type and line color
		self.pick_type = "P"
		self.pick_line_color = "red"

		# Initialise pick tracker
		self.pick_lines = {}

		# Initialise component rotation tracker
		self.components = "ZNE"

		# Initialise tracker for spectrogram window
		self.show_spectrogram = False

		# Initialise trace rejection tracker
		self.trace_removed = False

		# Initialise toggle trackers
		self.ctrl  = False
		self.shift = False

		# If picking all receivers for a given source
		if receiver == None and source != None:
			self.recs = True
		else:
			self.recs = False

		# If picking all sources at a given receiver
		if source == None and receiver != None:
			self.srcs = True
		else:
			self.srcs = False

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/trace_window.ui', self)
		self.uiToggleSpectrogramStacked.setCurrentIndex(0)

		# If just picking a single source at a single receiver
		if not self.srcs and not self.recs:
			# Disable the next trace function
			self.uiNextTraceButton.setEnabled(False)
			self.uiLastTraceButton.setEnabled(False)

		# If picking all sources at a given receiver
		elif self.srcs:
			# Load the sources
			sources = glob.glob('{}/data/{}/*.z'.format(self.catalogue.catalogue_path, receiver.upper()))
			self.sources = []
			for source in sources:
				head, tail  = os.path.split(source)
				self.sources.append(tail.split(".")[1])

			self.counter = 0

			try:
				self.source = psm.Source(self.sources[self.counter])
			except IndexError:
				qt.QMessageBox.about(self, "Error!", "There are no sources with recorded arrivals at this receiver.")
				return

		# If picking all receivers for a given source
		elif self.recs:
			# Load receivers
			receivers = glob.glob("{}/data/*/source.{}.*.z".format(self.catalogue.catalogue_path, self.source))
			self.receivers = []
			for receiver in receivers:
				head, tail = os.path.split(receiver)
				self.receivers.append(tail.split(".")[2])

			self.counter = 0

			try:
				self.receiver = self.receivers[self.counter]
			except IndexError:
				qt.QMessageBox.about(self, "Error!", "There are no receivers with recorded arrivals for this source.")
				return

			# Open up Wadati plot
			self.wadatiWindow = WadatiWindow(self)

		self._updateReceiverInformation()
		self._updateSourceInformation()

		# Initialise trace axes
		self.uiComponent1Mpl.canvas._tracePlot()
		self.uiComponent2Mpl.canvas._tracePlot()
		self.uiComponent3Mpl.canvas._tracePlot()

		# Initialise spectrogram axis
		self.uiSpectrogramPlot.canvas._tracePlot()

		if self.default_filt != None:
			self._populateFilter()

		self.plotTraces()

		self.connect()

		self.setWindowTitle('SeisPick - Trace picking window')
		self.setWindowIcon(QtGui.QIcon("misc/icon.png"))
		self.show()
	
	# ------------------------------

	# ------------------------------
	# Connection and event functions
	# ------------------------------
	
	def connect(self):
		# File menu connections
		self.uiSaveAction.triggered.connect(self.saveTrace)

		# Settings menu connections
		self.uiDefaultFilterAction.triggered.connect(self.defaultFilter)
		self.uiPPickAction.triggered.connect(lambda: self.updatePick("P"))
		self.uiSPickAction.triggered.connect(lambda: self.updatePick("S"))
		self.uiCustomPickAction.triggered.connect(lambda: self.updatePick("C"))

		# Connect to filter buttons
		self.cidafilter = self.uiApplyFilterButton.clicked.connect(self.applyFilter)
		self.cidrfilter = self.uiRemoveFilterButton.clicked.connect(self.removeFilter)

		# Connect to plot option buttons
		self.cidrotatecomps = self.uiRotateCompButton.clicked.connect(self.rotateComponents)
		self.cidspectrogram = self.uiViewSpectrogramButton.clicked.connect(self.toggleSpectrogram)
		self.cidnexttrace = self.uiNextTraceButton.clicked.connect(self.nextTrace)
		self.cidlasttrace = self.uiLastTraceButton.clicked.connect(self.previousTrace)
		self.cidrjcttrace = self.uiRejectTraceButton.clicked.connect(self.rejectTrace)
		self.cidlastzoom  = self.uiToggleLimitsButton.clicked.connect(self.toggleLims)
		self.cidresetplot = self.uiResetPlotButton.clicked.connect(self.resetPlot)
	
	def plotconnect(self):
		# Connect each canvas to track motion - this will create a 
		self.c1_moveid = self.uiComponent1Mpl.canvas.mpl_connect('motion_notify_event', self._onMove)
		self.c2_moveid = self.uiComponent2Mpl.canvas.mpl_connect('motion_notify_event', self._onMove)
		self.c3_moveid = self.uiComponent3Mpl.canvas.mpl_connect('motion_notify_event', self._onMove)

		# Connect to clicks
		self.c1_clickid = self.uiComponent1Mpl.canvas.mpl_connect('button_press_event', self._onClick)
		self.c2_clickid = self.uiComponent2Mpl.canvas.mpl_connect('button_press_event', self._onClick)
		self.c3_clickid = self.uiComponent3Mpl.canvas.mpl_connect('button_press_event', self._onClick)

		# Connect to releases
		self.c1_releaseid = self.uiComponent1Mpl.canvas.mpl_connect('button_release_event', self._onRelease)
		self.c2_releaseid = self.uiComponent2Mpl.canvas.mpl_connect('button_release_event', self._onRelease)
		self.c3_releaseid = self.uiComponent3Mpl.canvas.mpl_connect('button_release_event', self._onRelease)

	def plotdisconnect(self):
		# Disconnect each canvas to track motion - this will create a 
		self.uiComponent1Mpl.canvas.mpl_disconnect(self.c1_moveid)
		self.uiComponent2Mpl.canvas.mpl_disconnect(self.c2_moveid)
		self.uiComponent3Mpl.canvas.mpl_disconnect(self.c3_moveid)

		# Disconnect to clicks
		self.uiComponent1Mpl.canvas.mpl_disconnect(self.c1_clickid)
		self.uiComponent2Mpl.canvas.mpl_disconnect(self.c2_clickid)
		self.uiComponent3Mpl.canvas.mpl_disconnect(self.c3_clickid)

		# Disconnect to releases
		self.uiComponent1Mpl.canvas.mpl_disconnect(self.c1_releaseid)
		self.uiComponent2Mpl.canvas.mpl_disconnect(self.c2_releaseid)
		self.uiComponent3Mpl.canvas.mpl_disconnect(self.c3_releaseid)
	
	def saveTrace(self):
		if self.src != None:
			print("Saving source...")
			self.src.saveData()
		else:
			return

	def keyPressEvent(self, event):
		if event.key() == Qt.Key_Control:
			self.ctrl = True

		if event.key() == Qt.Key_Shift:
			self.shift = True

		if event.key() == Qt.Key_U or event.key() == Qt.Key_D:
			if not self.pick_line:
				return
			else:
				pick_polarity = chr(event.key())
				self.uiPolarityDisplay.setText(pick_polarity)
				self.src.addPick(info_type="polarity", 
								 value=pick_polarity, 
								 pick_type=self.pick_type)

	def keyReleaseEvent(self, event):
		# Toggle Ctrl Modifier off
		if event.key() == Qt.Key_Control:
			self.ctrl = False

		# Toggle Shift Modifier off
		if event.key() == Qt.Key_Shift:
			self.shift = False

	def _onMove(self, event):
		# Temp variables for accessing canvases
		c1_canvas = self.uiComponent1Mpl.canvas
		c2_canvas = self.uiComponent2Mpl.canvas
		c3_canvas = self.uiComponent3Mpl.canvas

		if event.inaxes is c1_canvas.ax or c2_canvas.ax or c3_canvas.ax:
			if self.shift:
				if not self._zoom_click:
					return
				if self._trace_drag_lock is not self:
					return

				xpress, ypress = self.trace_click

				xmove, ymove = event.xdata, event.ydata

				dx = xmove - xpress
				dy = ymove - ypress

			else:
				# Grab the x position of the pointer
				x = event.xdata

				# Set the x value of the cursor to the current position
				self.c1_cursor.set_data([x, x], self.c1_y)
				self.c2_cursor.set_data([x, x], self.c2_y)
				self.c3_cursor.set_data([x, x], self.c3_y)

			# restore the background region
			c1_canvas.restore_region(self.c1_bg)
			c2_canvas.restore_region(self.c2_bg)
			c3_canvas.restore_region(self.c3_bg)

			if self.shift:
				# Set rectangle values
				c1_rect = Rectangle((xpress, ypress), dx, dy,
										  edgecolor='red', fill=False)
				c2_rect = Rectangle((xpress, ypress), dx, dy,
										  edgecolor='red', fill=False)
				c3_rect = Rectangle((xpress, ypress), dx, dy,
										  edgecolor='red', fill=False)
				# Add the rectangle patch to each canvas
				c1_canvas.ax.add_patch(c1_rect)
				c2_canvas.ax.add_patch(c2_rect)
				c3_canvas.ax.add_patch(c3_rect)

				# Redraw just the rectangle
				c1_canvas.ax.draw_artist(c1_rect)
				c2_canvas.ax.draw_artist(c2_rect)
				c3_canvas.ax.draw_artist(c3_rect)

			else:
				# Redraw just the current cursor
				c1_canvas.ax.draw_artist(self.c1_cursor)
				c2_canvas.ax.draw_artist(self.c2_cursor)
				c3_canvas.ax.draw_artist(self.c3_cursor)

			# Plot any lines that are currently being stored
			self._replotLines()				

			# blit just the redrawn area
			c1_canvas.blit(c1_canvas.ax.bbox)
			c2_canvas.blit(c2_canvas.ax.bbox)
			c3_canvas.blit(c3_canvas.ax.bbox)

	def _onClick(self, event):
		# Check for shift modifier - this enables the draggable zoom
		if self.shift:
			if self._trace_drag_lock is not None:
				return

			self._zoom_click = True

			xpress, ypress = event.xdata, event.ydata
			self._trace_drag_lock = self

			# Grab the boundary values
			self.trace_click = xpress, ypress

		else:
			c1_canvas = self.uiComponent1Mpl.canvas
			c2_canvas = self.uiComponent2Mpl.canvas
			c3_canvas = self.uiComponent3Mpl.canvas

			dt = self.src.delta

			adjusted_xdata = round(event.xdata / dt) * dt

			# Left-clicking handles the window start time
			if event.button == 1 and not self.ctrl:
				# Set the window start line to be redrawn on move
				self.wbeg_vline = True

				# Add the window beginning to the event stats
				self.src.addPick(info_type="window", value=adjusted_xdata, pick_type="wb")

				# Make a vertical line artist
				self.c1_wbeg = c1_canvas.ax.axvline(adjusted_xdata, linewidth=1, color="green", animated=True)
				self.c2_wbeg = c2_canvas.ax.axvline(adjusted_xdata, linewidth=1, color="green", animated=True)
				self.c3_wbeg = c3_canvas.ax.axvline(adjusted_xdata, linewidth=1, color="green", animated=True)

				# Restore the background region
				c1_canvas.restore_region(self.c1_bg)
				c2_canvas.restore_region(self.c2_bg)
				c3_canvas.restore_region(self.c3_bg)

				# Plot any lines that are currently being stored
				self._replotLines()

				# blit the redrawn area
				c1_canvas.blit(c1_canvas.ax.bbox)
				c2_canvas.blit(c2_canvas.ax.bbox)
				c3_canvas.blit(c3_canvas.ax.bbox)

			# Middle-clicking handles the arrival pick time
			if event.button == 2 or (event.button == 1 and self.ctrl):
				# Set the pick line to be redrawn on move
				self.pick_line = True

				# Add the pick to the event stats
				self.src.addPick(info_type="pick", value=adjusted_xdata, pick_type=self.pick_type)

				# Set pick time label
				pick_time = self.src.starttime + adjusted_xdata
				self.uiPickTimeDisplay.setText(pick_time.isoformat())
				self.uiPhaseDisplay.setText(self.pick_type)

				# Make a vertical line artist
				self.comp_1_pick = c1_canvas.ax.axvline(adjusted_xdata, linewidth=1, color=self.pick_line_color, animated=True)
				self.comp_2_pick = c2_canvas.ax.axvline(adjusted_xdata, linewidth=1, color=self.pick_line_color, animated=True)
				self.comp_3_pick = c3_canvas.ax.axvline(adjusted_xdata, linewidth=1, color=self.pick_line_color, animated=True)

				self.pick_lines[self.pick_type] = [self.comp_1_pick, self.comp_2_pick, self.comp_3_pick]

				# Restore the background region
				c1_canvas.restore_region(self.c1_bg)
				c2_canvas.restore_region(self.c2_bg)
				c3_canvas.restore_region(self.c3_bg)

				# Plot any lines that are currently being stored
				self._replotLines()

				# blit the redrawn area
				c1_canvas.blit(c1_canvas.ax.bbox)
				c2_canvas.blit(c2_canvas.ax.bbox)
				c3_canvas.blit(c3_canvas.ax.bbox)

			# Right-clicking handles the window end time
			if event.button == 3:
				# Set the window end line to be redrawn on move
				self.w_end_line = True

				# Add the window ending to the event stats
				self.src.addPick(info_type="window", value=adjusted_xdata, pick_type="we")

				# Make a vertical line artist
				self.c1_wend = c1_canvas.ax.axvline(adjusted_xdata, linewidth=1, color="green", animated=True)
				self.c2_wend = c2_canvas.ax.axvline(adjusted_xdata, linewidth=1, color="green", animated=True)
				self.c3_wend = c3_canvas.ax.axvline(adjusted_xdata, linewidth=1, color="green", animated=True)

				# Restore the background region
				c1_canvas.restore_region(self.c1_bg)
				c2_canvas.restore_region(self.c2_bg)
				c3_canvas.restore_region(self.c3_bg)

				# Plot any lines that are currently being stored
				self._replotLines()

				# blit the redrawn area
				c1_canvas.blit(c1_canvas.ax.bbox)
				c2_canvas.blit(c2_canvas.ax.bbox)
				c3_canvas.blit(c3_canvas.ax.bbox)

	def _onRelease(self, event):
		if self.shift:
		
			if self._trace_drag_lock is not self:
				return

			# Grab the x and y data positions of the click point
			xpress, ypress = self.trace_click

			# Reset trace click and lock variables
			self._zoom_click      = False
			self.trace_click      = None
			self._trace_drag_lock = None

			# Grab the x and y data positions of the release point
			xrelease, yrelease = event.xdata, event.ydata

			# Set limits
			xlims = (min(xpress, xrelease), max(xpress, xrelease))
			ylims = (min(ypress, yrelease), max(ypress, yrelease))
			self.previous_lims = self.lims
			self.lims = [xlims, ylims]

			self.plotTraces(replot=True)

	# ------------------------------

	# ---------------------------
	# Handlers for trace plotting
	# ---------------------------

	def plotTraces(self, replot=False):

		c1_canvas = self.uiComponent1Mpl.canvas
		c2_canvas = self.uiComponent2Mpl.canvas
		c3_canvas = self.uiComponent3Mpl.canvas

		# Clear the canvases
		c1_canvas.ax.clear()
		c2_canvas.ax.clear()
		c3_canvas.ax.clear()

		# Create background variables
		self.c1_bg = None
		self.c2_bg = None
		self.c3_bg = None

		# Create the lock and click variables
		self._zoom_click      = False
		self._trace_drag_lock = None
		self.trace_click      = None

		if not replot:
			# Add limits variable
			self.previous_lims = None
			self.lims = None

			# Set filter to default filter
			self.filt = self.default_filt

			# Make the window and pick lines not show on move
			self.wbeg_vline = False
			self.pick_line  = False
			self.w_end_line = False

			# Create an instance of the Source class
			self.src = psm.SourceReceiverPair(self.catalogue_path, self.source, self.receiver)

			# Look up any picks and make vertical lines
			for pick, data in self.src.picks.items():
				print(pick)
				pick_time = data["rtime"]
				if pick == "we" or pick == "wb":
					color = "green"
				else:
					color = "gray"

				c1_canvas.ax.axvline(pick_time, linewidth=1, color=color)
				c2_canvas.ax.axvline(pick_time, linewidth=1, color=color)
				c3_canvas.ax.axvline(pick_time, linewidth=1, color=color)

		else:
			self._replotLines()

		# Check if a filter has been specified, and apply it if so
		if not self.filt == None:
			self.src.filterObspy(self.filt)

		# Plot the traces
		self.src.plotTraces([c1_canvas.ax, c2_canvas.ax, c3_canvas.ax], lims=self.lims)

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
		# Grab the trace background when it is drawn
		self.c1_bg = event.canvas.copy_from_bbox(event.canvas.ax.bbox)

	def _c2DrawEvent(self, event):
		# Grab the trace background when it is drawn
		self.c2_bg = event.canvas.copy_from_bbox(event.canvas.ax.bbox)

	def _c3DrawEvent(self, event):
		# Grab the trace background when it is drawn
		self.c3_bg = event.canvas.copy_from_bbox(event.canvas.ax.bbox)

	def _replotLines(self):
		# For clarity
		c1_canvas = self.uiComponent1Mpl.canvas
		c2_canvas = self.uiComponent2Mpl.canvas
		c3_canvas = self.uiComponent3Mpl.canvas

		if self.wbeg_vline:
			# Redraw the window line
			c1_canvas.ax.draw_artist(self.c1_wbeg)
			c2_canvas.ax.draw_artist(self.c2_wbeg)
			c3_canvas.ax.draw_artist(self.c3_wbeg)

		if self.pick_line:
			# Redraw the pick lines
			for phase, pick_lines in self.pick_lines.items():
				c1_canvas.ax.draw_artist(pick_lines[0])
				c2_canvas.ax.draw_artist(pick_lines[1])
				c3_canvas.ax.draw_artist(pick_lines[2])

		if self.w_end_line:
			# Redraw the window line
			c1_canvas.ax.draw_artist(self.c1_wend)
			c2_canvas.ax.draw_artist(self.c2_wend)
			c3_canvas.ax.draw_artist(self.c3_wend)

	def _updateReceiverInformation(self):
		self.uiReceiverNameDisplay.setText(self.receiver.station)
		self.uiReceiverLonDisplay.setText(f"{self.receiver.longitude:.4f}")
		self.uiReceiverLatDisplay.setText(f"{self.receiver.latitude:.4f}")
		self.uiReceiverElevDisplay.setText(f"{self.receiver.elevation:.4f}")
		self.uiReceiverDepDisplay.setText(self.receiver.deployment.isoformat().split("T")[0])
		self.uiReceiverRetDisplay.setText(self.receiver.retrieval.isoformat().split("T")[0])

	def _updateSourceInformation(self):
		self.uiOriginDateDisplay.setText(self.source.otime.isoformat().split("T")[0])
		self.uiOriginTimeDisplay.setText(self.source.otime.isoformat().split("T")[1])
		self.uiSourceLonDisplay.setText(f"{self.source.longitude:.4f}")
		self.uiSourceLatDisplay.setText(f"{self.source.latitude:.4f}")
		self.uiSourceDepthDisplay.setText(f"{self.source.depth:.4f}")
		if (type(self.source.magnitude) == float) or (type(self.source.magnitude) == np.float64):
			self.uiSourceMagDisplay.setText(f"{self.source.magnitude:.2f}")
		else:
			self.uiSourceMagDisplay.setText(str(self.source.magnitude))
		self.uiSourceIDDisplay.setText(str(self.source.sourceid))

	# ---------------------------

	# ---------------
	# Filter handlers
	# ---------------

	def defaultFilter(self):
		# Open default filter argument as a child of parent class
		self.parent.defaultFilter()

	def applyFilter(self):
		# Read the default filter parameters
		try:
			filt_type = self.uiFilterTypeComboBox.currentText()
			no_poles  = int(self.uiPoleCountComboBox.currentText())
			zerophase = self.uiZeroPhaseCheckBox.isChecked()
			low_freq  = float(self.uiLowFreqInput.text())
			high_freq = float(self.uiHighFreqInput.text())

			self.filt = {'filt_type': filt_type,
						 'no_poles': no_poles,
						 'zerophase': zerophase,
						 'low_freq': low_freq,
						 'high_freq': high_freq
						}
		except ValueError:
			qt.QMessageBox.about(self, "Warning!", "You appear to have tried to use an incomplete/incorrect filter. Fill in all of the options and try again.")
			return

		self.plotdisconnect()
		self.plotTraces(replot=True)

	def removeFilter(self):
		self.plotdisconnect()
		self.filt = None
		self.src.removeFilter()
		self.plotTraces(replot=True)

	def _populateFilter(self):
		# Read in default filter type and options
		filter_type    = self.default_filt[0]
		filter_options = self.default_filt[1]

		self.uiFilterTypeComboBox.setCurrentText(filter_type)
		self.uiPoleCountComboBox.setCurrentText(str(filter_options["corners"]))

		if filter_options["zerophase"]:
			self.uiZeroPhaseCheckBox.setChecked(True)

		if filter_type == "bandpass":
			self.uiLowFreqInput.setText(str(filter_options["minfreq"]))
			self.uiHighFreqInput.setText(str(filter_options["maxfreq"]))

		if filter_type == "lowpass":
			self.uiLowFreqInput.setText(str(filter_options["freq"]))

		if filter_type == "highpass":
			self.uiHighFreqInput.setText(str(filter_options["freq"]))

	# ---------------

	# -------------------
	# Phase pick handler
	# -------------------

	def updatePick(self, pick_type):
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

			# Open custom phase dialogue
			self.customPickDialogue = CustomPickDialogue()

			if not self.customPickDialogue.exec_():
				return

		self.uiPolarityDisplay.setText("")
	
	# -------------------

	# ---------------------
	# Plot actions handlers
	# ---------------------
	def rotateComponents(self):
		if self.components == "ZNE":
			# Switch component tracker
			self.components = "ZRT"

			# Rotate components
			self.src.rotate(method="NE->RT")

		elif self.components == "ZRT":
			# Switch component tracker
			self.components = "ZNE"

			# Rotate components
			self.src.rotate(method="RT->NE")

		self.plotTraces(replot=True)

	def toggleSpectrogram(self):
		self.show_spectrogram = not self.show_spectrogram
		self.plotSpectrogram()

	def plotSpectrogram(self):
		if self.show_spectrogram:
			self.uiToggleSpectrogramStacked.setCurrentIndex(1)
			self.src.plotSpectrogram(self.uiSpectrogramPlot.canvas.ax)

		else:
			self.uiToggleSpectrogramStacked.setCurrentIndex(0)

	def nextTrace(self):
		# If previous trace was removed, do not advance counter
		if not self.trace_removed:
			# Advance counter
			self.counter += 1

		self.updateTrace()

	def previousTrace(self):
		# Reverse counter
		self.counter += -1

		self.updateTrace()

	def updateTrace(self):
		# Check if P and S have been picked (and window exists)
		if "P_manual" in self.src.picks.keys() and "S_manual" in self.src.picks.keys() and self.recs:
			ptravel = self.src.starttime - self.src.otime
			ptravel += self.src.picks["P_manual"]["rtime"]
			stravel = self.src.starttime - self.src.otime
			stravel += self.src.picks["S_manual"]["rtime"]
			self.wadatiWindow.addPick(ptravel, stravel, self.receiver)

		# Save any picks
		self.saveTrace()

		# Reset to P pick
		self.uiPRadio.setChecked(True)

		# Reset trace removal tracker
		self.trace_removed = False

		# Reset pick tracker
		self.pick_lines = {}
		self.pick_time = ""
		self.pick_phase = ""

		# Set all the labels
		self.uiPolarityDisplay.setText("")
		self.uiPickTimeDisplay.setText(self.pick_time)
		self.uiPhaseDisplay.setText(self.pick_phase)

		if self.srcs:
			# Check not trying to go below 0th position
			if self.counter == -1:
				qt.QMessageBox.about(self, "Warning!", "Looped round to last source.")
				self.counter = len(sources) - 1

			# Check not on last source
			elif self.counter == len(self.sources):
				qt.QMessageBox.about(self, "Warning!", "Looped round to first source.")
				self.counter = 0

			# Change current source
			self.source = self.sources[self.counter]

			# Update source information
			self._updateSourceInformation(self.source)

		if self.recs:
			# Check not trying to go below 0th position
			if self.counter == -1:
				qt.QMessageBox.about(self, "Warning!", "Looped round to last receiver.")
				self.counter = len(receivers) - 1

			# Check not on last receiver
			elif self.counter == len(self.receivers):
				qt.QMessageBox.about(self, "Warning!", "Looped round to first receiver.")
				self.counter = 0

			# Change current receiver
			self.receiver = self.receivers[self.counter]

			# Update receiver information
			self._updateReceiverInformation(self.receiver)

		self.plotdisconnect()
		self.plotTraces()

	def rejectTrace(self):
		if self.src != None:
			print("Deleting trace...")
			traces = glob.glob(self.src.file_path)
			for trace in traces:
				os.remove(trace)

			# Pass information back to parent and remove arrival from available arrivals?

			# Delete source or receiver from sources/receivers list
			if self.srcs:
				del self.sources[self.counter]

			if self.recs:
				del self.receivers[self.counter]

			self.trace_removed = True

			self.nextTrace()
		else:
			return

	def toggleLims(self):
		# Overwrite limits
		self.lims, self.previous_lims = self.previous_lims, self.lims

		# Replot the traces
		self.plotTraces(replot=True)

	def resetPlot(self):
		self.plotdisconnect()
		self.lims = None
		self.filt = self.default_filter
		self.plotTraces()

	# ---------------------


if __name__ == "__main__":

	app = qt.QApplication(sys.argv)
	seispick = SeisPick()
	sys.exit(app.exec_())