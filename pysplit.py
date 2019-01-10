#!/usr/bin/env python3

"""
The PySplit GUI class system - this class handles the connection between
the backend (Event and Catalogue class systems) and the user interface.

Methods
 

Attributes


Author: Hemmelig
"""

import numpy as np
import pandas as pd
from mpl_toolkits.basemap import Basemap
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import sys
from PyQt5.QtCore import Qt, QDate, QDateTime
import PyQt5.QtWidgets as qt
import PyQt5.QtGui as QtGui
from obspy import read
from PyQt5 import uic
from mpl_toolkits.basemap import Basemap
from string import Template
import glob
import os
from obspy import UTCDateTime

import catalogue as cat
import event as evt

class PySplit(qt.QMainWindow):
	"""
	The PySplit class - base class for generation of, and interaction with, the PySplit GUI.
	"""

	def __init__(self):
		super(PySplit, self).__init__()

		self.initUI()

	def initUI(self):

		# Load the ui file
		uic.loadUi('ui_files/main.ui', self)

		# Initialise toggle trackers
		self.ctrl_toggle  = False
		self.shift_toggle = False

		# Connect to actions and interactive widgets
		self.connect()

		# Hide the map widget initially
		self.mpl.hide()
		self.catalogue_forms.setCurrentIndex(2)

		# Set initial filter to None
		self.filt = None

		# Specify input options and disable them on load
		self.input_options = [self.startDate_input, self.endDate_input, self.minMag_input, self.centreLon_input, self.centreLat_input, self.minRad_input, self.maxRad_input, self.tele_load_button]

		self.setWindowTitle('PySplit - Shear Wave Splitting software')
		self.setWindowIcon(QtGui.QIcon("misc/icon.png"))
		self.show()

	def connect(self):
		# File menu connections
		self.actionNew_catalogue.triggered.connect(self.newCatalogue)
		self.actionLoad_catalogue.triggered.connect(self.loadCatalogue)

		# Settings menu connections
		self.actionDefault_filter.triggered.connect(self.defaultFilter)

		# Map connections
		self.mpl.canvas.mpl_connect('pick_event', self._onPick)
		self.mpl.canvas.mpl_connect('button_press_event', self._onMapClick)
		self.mpl.canvas.mpl_connect('motion_notify_event', self._onMapMove)
		self.mpl.canvas.mpl_connect('button_release_event', self._onMapRelease)

		# Station list connection
		self.station_list.doubleClicked.connect(self.stationSelect)

		# Event list connection
		self.events_list.doubleClicked.connect(self.eventSelect)

		# Catalogue generation button connection
		self.tele_load_button.clicked.connect(self.generateTeleseismicCatalogue)

		# Catalogue load button connection
		self.load_events_button.clicked.connect(self.load_events)
		self.load_stats_button.clicked.connect(self.load_stations)
		self.load_arrivals_button.clicked.connect(self.load_arrivals)
		self.load_waveforms_button.clicked.connect(self.load_waveforms)

		# Plot option tickboxes connection
		self.picked_plot_tickbox.stateChanged.connect(self.pickStateChanged)
		self.unpicked_plot_tickbox.stateChanged.connect(self.unpickedStateChanged)
		self.station_plot_tickbox.stateChanged.connect(self.stationStateChanged)

		self.button_plotMap.clicked.connect(lambda: self.plotCatalogueMap(replot=True))
		self.button_resetMap.clicked.connect(self.plotCatalogueMap)

	def keyPressEvent(self, event):
		# Toggle Ctrl Modifier on
		if event.key() == Qt.Key_Control:
			self.ctrl_toggle = True

		# Toggle Shift Modifier on
		if event.key() == Qt.Key_Shift:
			self.shift_toggle = True

	def keyReleaseEvent(self, event):
		# Toggle Ctrl Modifier off
		if event.key() == Qt.Key_Control:
			self.ctrl_toggle = False

		# Toggle Shift Modifier off
		if event.key() == Qt.Key_Shift:
			self.shift_toggle = False

	def unpickedStateChanged(self, int):
		if self.unpicked_plot_tickbox.isChecked():
			pass
			#self.catalogue.station_points.set_visible(True)
		else:
			pass
			#self.catalogue.station_points.set_visible(False)

	def stationStateChanged(self, int):
		if self.station_plot_tickbox.isChecked():
			pass
			#self.catalogue.station_points.set_visible(True)
		else:
			pass
			#self.catalogue.station_points.set_visible(False)

	def pickStateChanged(self, int):
		if self.picked_plot_tickbox.isChecked():
			pass
			#self.catalogue.pick_points.set_visible(True)
		else:
			pass
			#self.catalogue.pick_points.set_visible(False)

	def _onPick(self, event):
	    artist = event.artist
	    xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
	    xy = artist.get_offsets()
	    label = artist.get_label()
	    ind = event.ind
	    lon, lat = self.catalogue.m(xy[ind[0]][0], xy[ind[0]][1], inverse=True)

	    # If station
	    if "STAT" in label:
	    	stat = label.split(": ")[1]
	    	self._updateStationInformation(stat)

	    	# If in station pick mode, open a Picking Window
	    	if self.actionPick_stations.isChecked():
	    		self.pickWindow = PickingWindow(self.catalogue, self.catalogue_name, filt=self.filt, station=stat)

	    # If event
	    if "EVENT" in label:
	    	evt = label.split(": ")[1]
	    	self._updateEventInformation(evt)

	    	# If in event pick mode, open a Picking Window
	    	if self.actionPick_events.isChecked():
		    	self.pickWindow = PickingWindow(self.catalogue, self.catalogue_name, filt=self.filt, event=evt)

	def defaultFilter(self):
		self.defaultFilterDialogue = DefaultFilterDialogue()

		if self.defaultFilterDialogue.exec_():
			# Read the default filter parameters
			try:
				filt_type = self.defaultFilterDialogue.input_filtType.currentText()
				no_poles  = int(self.defaultFilterDialogue.input_noPoles.currentText())
				zerophase = self.defaultFilterDialogue.input_zerophase.isChecked()
				low_freq  = float(self.defaultFilterDialogue.input_lowFreq.text())
				high_freq = float(self.defaultFilterDialogue.input_highFreq.text())

				self.filt = {'filt_type': filt_type,
							 'no_poles': no_poles,
							 'zerophase': zerophase,
							 'low_freq': low_freq,
							 'high_freq': high_freq
							}
			except:
				qt.QMessageBox.about(self, "Error!", "Try again.")
		else:
			return

	def newCatalogue(self):
		self.newCatalogueDialogue = NewCatalogueDialogue(self)

		if not self.newCatalogueDialogue.exec_():
			return

		if self.catalogue_type == "teleseismic":
			# Set catalogue type and data source labels
			self.data_source = "IRIS"
			self.label_catTypeDisp.setText("Teleseismic")

			# Display teleseismic catalogue parameter form page
			self.catalogue_forms.setCurrentIndex(1)
			self.plot_options.setCurrentIndex(1)

			# Create an instance of the Teleseismic Catalogue class
			self.catalogue = cat.TeleseismicCatalogue("{}/{}".format(self.catalogue_path, self.catalogue_name), self.archive_path, self.receiver_file)

		elif self.catalogue_type == "local":
			# Set catalogue type label
			self.label_catTypeDisp.setText("Local")

			# Display local catalogue progress page
			self.catalogue_forms.setCurrentIndex(0)
			self.plot_options.setCurrentIndex(0)

			# Get input file
			self.localInputDialogue = LocalInputDialogue(self)
			if not self.localInputDialogue.exec_():
				return

			# Create an instance of the Local Catalogue class
			self.catalogue = cat.LocalCatalogue("{}/{}".format(self.catalogue_path, self.catalogue_name), self.archive_path, self.receiver_file)

			# Create the catalogue metafile - this contains general information about the catalogue
			# and is read in when an existing catalogue is loaded
			self._generateCatalogueMetafile()
			self.generate_catalogue()

			# Plot the map of the catalogue
			self.plotCatalogueMap()

		# Set the catalogue name label
		self.label_catNameDisp.setText(self.catalogue_name)

		# Populate the station list
		self._populateStationList()

		# Enable the station and event selection lists
		self.station_list.setEnabled(True)
		self.events_list.setEnabled(True)

	def _populateStationList(self):
		# Populate the station list
		model = QtGui.QStandardItemModel(self.station_list)
		for receiver in self.catalogue.receiver_df.stat.values:
			item = QtGui.QStandardItem(receiver)
			model.appendRow(item)
		self.station_list.setModel(model)
		self.station_list.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)

	def _populateCatalogueInformation(self):
		self.label_startDate2.setText(self.sdate)
		self.label_endDate2.setText(self.edate)
		self.label_cenLon.setText(self.clon)
		self.label_cenLat.setText(self.clat)
		self.label_minRad2.setText(self.minrad)
		self.label_maxRad2.setText(self.maxrad)
		self.label_minMag2.setText(self.minmag)

		# Calculate number of events and update label
		no_events = str(len(self.catalogue.source_df.index))
		self.label_noEvents.setText(no_events)

	def plotCatalogueMap(self, replot=False):
		# Transmit a message to the status bar
		self.statusbar.showMessage("Plotting catalogue map...")

		# Create background, lock and click variables
		self.map_background = None
		self._map_drag_lock = None
		self.map_click      = None

		# Show the matplotlib widget and clear it
		self.mpl.show()
		self.mpl.canvas.ax.clear()

		# If replotting, grab the new parameters
		if replot:
			self.catalogue.plot_geographic(self.mpl.canvas.ax,
										   lon0=float(self.input_minLon.text()),
										   lon1=float(self.input_maxLon.text()),
										   lat0=float(self.input_minLat.text()),
										   lat1=float(self.input_maxLat.text()))

			if self.station_plot_tickbox.isChecked():
				self.catalogue.plot_stations(self.mpl.canvas.ax)

		else:
			# Plot the events
			self.catalogue.plot_geographic(self.mpl.canvas.ax)

			# Plot the stations
			if self.station_plot_tickbox.isChecked():
				self.catalogue.plot_stations(self.mpl.canvas.ax)

			# Grab the latitude and longitude of the map and send to input options
			if self.catalogue_type == "local":
				self.input_minLon.setText(str(f"{self.catalogue.lon0:.5f}"))
				self.input_maxLon.setText(str(f"{self.catalogue.lon1:.5f}"))
				self.input_minLat.setText(str(f"{self.catalogue.lat0:.5f}"))
				self.input_maxLat.setText(str(f"{self.catalogue.lat1:.5f}"))

		# Connect to the map to grab background once Qt has done resizing
		self.mpl.canvas.mpl_connect('draw_event', self._mapDrawEvent)

		# Draw the canvas once Qt has done all resizing
		self.mpl.canvas.draw_idle()

		# Transmit a message to the status bar
		self.statusbar.showMessage("Catalogue load complete.")

	def _mapDrawEvent(self, evt):
		self.map_background = evt.canvas.copy_from_bbox(evt.canvas.ax.bbox)

	def stationSelect(self, index):
		# Parse station selected
		self.station = self.station_list.model().data(index)

		# Load station information and print to display
		self._updateStationInformation(self.station)

		# Filter and load the corresponding events
		station_path = "{}/{}/data/{}".format(self.catalogue_path, self.catalogue_name, self.station.upper())

		events = glob.glob('{}/*.z'.format(station_path))

		# Populate the event list
		model = QtGui.QStandardItemModel(self.events_list)
		for event in events:
			head, tail = os.path.split(event)
			item = QtGui.QStandardItem("Event {}".format(tail.split(".")[1]))
			model.appendRow(item)
		self.events_list.setModel(model)
		self.events_list.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)

	def eventSelect(self, index):
		# Parse event selected
		self.event = self.events_list.model().data(index).split(" ")[1]

		# Load event information and print to display
		self._updateEventInformation(self.event)

		# Open up picking window
		self.pickWindow = PickingWindow(self.catalogue, self.catalogue_name, filt=self.filt, station=self.station, event=self.event)

	def _updateStationInformation(self, station):
		station_info = self.catalogue.receiver_df.query('stat == @station')
		self.label_statName.setText(station)
		self.label_statLon.setText(f"{station_info.lon.values[0]:.4f}")
		self.label_statLat.setText(f"{station_info.lat.values[0]:.4f}")
		self.label_statElev.setText(f"{station_info.dep.values[0]:.4f}")
		self.label_statDep.setText(station_info.st_dep.values[0].isoformat().split("T")[0])
		self.label_statRet.setText(station_info.et_dep.values[0].isoformat().split("T")[0])

	def _updateEventInformation(self, event):
		event_info = self.catalogue.source_df.query('sourceid == @event')
		self.label_eventOdate.setText(event_info.otime.values[0].isoformat().split("T")[0])
		self.label_eventOtime.setText(event_info.otime.values[0].isoformat().split("T")[1])
		self.label_eventLon.setText(f"{event_info.evlon.values[0]:.4f}")
		self.label_eventLat.setText(f"{event_info.evlat.values[0]:.4f}")
		self.label_eventDepth.setText(f"{event_info.evdep.values[0]:.4f}")
		if (type(event_info.evmag.values[0]) == float) or (type(event_info.evmag.values[0]) == np.float64):
			self.label_eventMag.setText(f"{event_info.evmag.values[0]:.2f}")
		else:
			self.label_eventMag.setText(str(event_info.evmag.values[0]))
		self.label_eventID.setText(str(event_info.sourceid.values[0]))

	def loadCatalogue(self):
		# Open file dialogue and choose catalogue directory
		self.catalogue_path = qt.QFileDialog.getExistingDirectory(self, 'Choose catalogue directory')

		# Test validity of the catalogue path
		if (os.path.exists("{}/data".format(self.catalogue_path))) and (os.path.exists("{}/metafiles".format(self.catalogue_path))):
			pass
		else:
			qt.QMessageBox.about(self, "Error!", "The directory you have chosen is invalid, please try again.")
			return

		# Parse the catalogue metafile
		try:
			self._parse_catalogue_metafile()
		except FileNotFoundError:
			qt.QMessageBox.about(self, "Error!", "Unable to find file containing the metadata for this catalogue. Please double-check it exists.")
			return

		# Set the catalogue name label
		self.label_catNameDisp.setText(self.catalogue_name)
		self.label_catTypeDisp.setText(self.catalogue_type)

		# Create an instance of Catalogue and load in the requisite information
		if self.catalogue_type == "local":
			# Create an instance of Local Catalogue
			self.statusbar.showMessage("Creating an instance of LocalCatalogue...")
			self.catalogue = cat.LocalCatalogue("{}/{}".format(self.catalogue_path, self.catalogue_name), self.archive_path, self.receiver_file)

			# Show the status page catalogues
			self.catalogue_forms.setCurrentIndex(0)
			self.plot_options.setCurrentIndex(0)

		if self.catalogue_type == "teleseismic":
			# Create an instance of Teleseismic Catalogue
			self.statusbar.showMessage("Creating an instance of TeleseismicCatalogue...")
			self.catalogue = cat.TeleseismicCatalogue("{}/{}".format(self.catalogue_path, self.catalogue_name), self.archive_path, self.receiver_file)

			# Show the catalogue info page
			self.catalogue_forms.setCurrentIndex(0)
			self.plot_options.setCurrentIndex(1)

		# Load the catalogue
		self.statusbar.showMessage("Loading catalogue...")
		self.load_events()
		self.load_stations()

		if self.catalogue_type == "teleseismic":
			# Display information about the catalogue in the catalogue info box
			self._populateCatalogueInformation()

		# Populate the station list
		self._populateStationList()

		# Enable the station and event selection lists
		self.station_list.setEnabled(True)
		self.events_list.setEnabled(True)

		# Plot the map of the catalogue
		self.plotCatalogueMap()

	def generateTeleseismicCatalogue(self):
		self.generate_catalogue()

	def load_events(self):
		sources = self.catalogue.load_sources()

		if sources:
			self.label_eventStat.setText("Events generated.")
			self.label_eventStat.setStyleSheet("color: rgb(0, 255, 0)")
			self.load_events_button.setEnabled(False)
		else:
			qt.QMessageBox.about(self, "Error!", "Error loading event file.")
			return

	def load_stations(self):
		stations = self.catalogue.load_receivers()

		if stations:
			self.label_statStat.setText("Stations generated.")
			self.label_statStat.setStyleSheet("color: rgb(0, 255, 0)")
			self.load_stats_button.setEnabled(False)
		else:
			qt.QMessageBox.about(self, "Error!", "Error loading station file.")
			return

	def load_arrivals(self):
		arrivals = self.catalogue.load_arrivals()
		print("1")

		if arrivals:
			self.label_arrivalStat.setText("Arrivals generated.")
			self.label_arrivalStat.setStyleSheet("color: rgb(0, 255, 0)")
			self.load_arrivals_button.setEnabled(False)
			self.load_waveforms_button.setEnabled(True)
		else:
			self.statusbar.showMessage("The arrivals file does not exist - attempting to generate one...")
			if self.catalogue_type == "local":
				# Load the arrivals from the local input file provided
				self.catalogue.get_arrivals(input_file=self.local_input, input_type=self.data_source)

			elif self.catalogue_type == "teleseismic":
				# Open pop-up to get the phases to be collected
				self.telePhaseDialogue = TelePhaseDialogue()

				if self.telePhaseDialogue.exec_():
					self.input_phases = self.telePhaseDialogue.phases
					if not self.input_phases:
						qt.QMessageBox.about(self, "Error!", "No phases selected - please try again!")
						return
				else:
					return

				self.catalogue.get_arrivals(phases=self.input_phases)
				self._populateStationList()

			self.load_arrivals()

	def load_waveforms(self):
		waveforms = self.catalogue.load_waveforms()

		if waveforms:
			self.label_waveformStat.setText("Waveforms downloaded.")
			self.label_waveformStat.setStyleSheet("color: rgb(0, 255, 0)")
			self.load_waveforms_button.setEnabled(False)
		else:
			self.catalogue.get_waveforms()
			#qt.QMessageBox.about(self, "Error!", "Some error in there somewhere yo")
			return

		self.label_waveformStat.setText("Waveforms generated.")
		self.label_waveformStat.setStyleSheet("color: rgb(0, 255, 0)")
		self.load_waveforms_button.setEnabled(False)

	def generate_catalogue(self):

		# Get catalogue type
		if self.catalogue_type == "local":
			self.catalogue.generate_catalogue(input_file=self.local_input, input_type=self.data_source)

		if self.catalogue_type == "teleseismic":
			self.catalogue.generate_catalogue(starttime=self.startDate_input.date().toString(Qt.ISODate), 
											  endtime=self.endDate_input.date().toString(Qt.ISODate),
											  minmag=float(self.minMag_input.text()),
											  lon=float(self.centreLon_input.text()),
											  lat=float(self.centreLat_input.text()),
											  minrad=float(self.minRad_input.text()),
											  maxrad=float(self.maxRad_input.text()))

			# Create the catalogue metafile - this contains general information about the catalogue
			# and is read in when an existing catalogue is loaded
			self._generateCatalogueMetafile()

			self.catalogue_forms.setCurrentIndex(0)

			self._populateCatalogueInformation()

			# Plot the map of the catalogue
			self.plotCatalogueMap()

	def _generateCatalogueMetafile(self):
		# Read in the template file
		filein = open('metafile_template.txt')
		src = Template(filein.read())

		self.sdate = self.startDate_input.date().toString(Qt.ISODate)
		self.edate = self.endDate_input.date().toString(Qt.ISODate)

		# Set the general metadata parameters
		d_cat = {'name': self.catalogue_name,
			 'type': self.catalogue_type,
			 'cpath': self.catalogue_path,
			 'dsource': self.data_source,
			 'cdate': QDateTime.currentDateTime().toUTC().toString(Qt.ISODate),
			 'apath': self.archive_path,
			 'aformat': self.archive_format,
			 'rfile': self.receiver_file,
			 'sdate': self.sdate,
			 'edate': self.edate}

		# Set the local catalogue parameters
		if self.catalogue_type == "local":
			d_loc = {'lfile': self.local_input}

			# Merge the two catalogues
			d = {**d_cat, **d_loc}

		# Set the teleseismic catalogue parameters
		if self.catalogue_type == "teleseismic":
			self.minmag = self.minMag_input.text()
			self.clon   = self.centreLon_input.text()
			self.clat   = self.centreLat_input.text()
			self.minrad = self.minRad_input.text()
			self.maxrad = self.maxRad_input.text()

			d_tel = {'minmag': self.minmag,
					 'clon': self.clon,
		 			 'clat': self.clat,
					 'minrad': self.minrad,
		 			 'maxrad': self.maxrad}

		 	# Merge the two catalogues
			d = {**d_cat, **d_tel}

		# Write out the dictionaries of the parameters - safe substitute ignores missing template variables
		# and just writes a blank
		output = src.safe_substitute(d)

		# Write out the catalogue metafile
		with open('{}/{}/metafiles/catalogue_metafile.txt'.format(self.catalogue_path, self.catalogue_name), 'w') as o:
			o.write(output)

	def _parse_catalogue_metafile(self):
		# Update status bar
		self.statusbar.showMessage('Parsing catalogue metafile...')
		# Empty list for collecting parameters
		params = []

		# Read in the parameters from the file
		for line in open("{}/metafiles/catalogue_metafile.txt".format(self.catalogue_path), 'r'):
			if "?" in line:
				line = line.split(" ? ")
				if line[1].endswith("\n"):
					param = line[1][:-1]
				params.append(param)

		# Parse general catalogue parameters
		self.catalogue_name = params[0]
		self.catalogue_type = params[1]
		self.catalogue_path = params[2]
		self.data_source    = params[3]
		self.cdate          = params[4]
		self.archive_path   = params[5]
		self.archive_format = params[6]
		self.receiver_file  = params[7]
		self.sdate          = params[8]
		self.edate          = params[9]

		# Parse local parameters
		if self.catalogue_type == "local":
			self.local_input = params[10]

		# Parse teleseismic parameters
		if self.catalogue_type == "teleseismic":
			self.minmag = params[11]
			self.clon   = params[12]
			self.clat   = params[13]
			self.minrad = params[14]
			self.maxrad = params[15]

	def _onMapClick(self, event):
		if not self.shift_toggle:
			return
		if self.catalogue_type == "teleseismic":
			return
		# See if mouse over map
		if event.inaxes != self.mpl.canvas.ax:
			return
		if self._map_drag_lock is not None:
			return

		xpress, ypress = event.xdata, event.ydata
		self._map_drag_lock = self

		# Grab the lon lat values 
		self.map_click = xpress, ypress

	def _onMapMove(self, event):
		if not self.shift_toggle:
			return
		if self.catalogue_type == "teleseismic":
			return

		if self._map_drag_lock is not self:
			return

		if event.inaxes != self.mpl.canvas.ax:
			return

		xpress, ypress = self.map_click

		xmove, ymove = event.xdata, event.ydata

		dx = xmove - xpress
		dy = ymove - ypress

		# Draw background from pixel buffer
		self.mpl.canvas.restore_region(self.map_background)

		# Set rectangle values
		map_rectangle = Rectangle((xpress, ypress), dx, dy,
								  edgecolor='red', fill=False)
		self.mpl.canvas.ax.add_patch(map_rectangle)
		self.mpl.canvas.ax.draw_artist(map_rectangle)

		# Blit the redrawn area
		self.mpl.canvas.blit(self.mpl.canvas.ax.bbox)

	def _onMapRelease(self, event):
		if not self.shift_toggle:
			return
		if self.catalogue_type == "teleseismic":
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
		self.mpl.canvas.restore_region(self.map_background)

		# Draw the final rectangle
		map_rectangle = Rectangle((xpress, ypress), dx, dy,
								  edgecolor='red', fill=False)
		self.mpl.canvas.ax.add_patch(map_rectangle)
		self.mpl.canvas.ax.draw_artist(map_rectangle)

		# Blit the redrawn area
		self.mpl.canvas.blit(self.mpl.canvas.ax.bbox)

		# Convert the x and y positions to lon/lat
		lonpress, latpress     = self.catalogue.m(xpress, ypress, inverse=True)
		lonrelease, latrelease = self.catalogue.m(xrelease, yrelease, inverse=True) 

		# Set the text values of the lon/lat input boxes
		self.input_minLon.setText(str(f"{min(lonpress, lonrelease):.2f}"))
		self.input_maxLon.setText(str(f"{max(lonpress, lonrelease):.2f}"))
		self.input_minLat.setText(str(f"{min(latpress, latrelease):.2f}"))
		self.input_maxLat.setText(str(f"{max(latpress, latrelease):.2f}"))

class NewCatalogueDialogue(qt.QDialog):

	def __init__(self, parent):
		super(NewCatalogueDialogue, self).__init__()

		self.parent = parent

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/new_catalogue_dialogue.ui', self)

		self.catDirButton.clicked.connect(self.browseCatalogue)
		self.arcDirButton.clicked.connect(self.browseArchive)
		self.receiverFileButton.clicked.connect(self.browseReceivers)
		self.buttonBox.accepted.connect(self.actionAccept)
		self.buttonBox.rejected.connect(self.actionReject)

		self.setWindowTitle('PySplit - New catalogue')
		self.show()

	def browseCatalogue(self):

		self.parent.label_eventOdate.setText("test")

		self.catalogue_path = qt.QFileDialog.getExistingDirectory(self, 'Choose catalogue directory')

		self.catBox.setText(self.catalogue_path)

	def browseArchive(self):

		self.archive_path = qt.QFileDialog.getExistingDirectory(self, 'Choose archive directory')

		self.arcBox.setText(self.archive_path)

	def browseReceivers(self):

		filename = qt.QFileDialog.getOpenFileName(self, 'Open file')

		self.receiver_file = filename[0]

		self.recBox.setText(self.receiver_file)

	def actionAccept(self):
		# If the inputs are accepted, set all of the parameters within the parent class
		# Test if any of the input fields have not been completed
		if (self.catBox.text() == "") or (self.catName.text() == "") or (self.recBox.text() == ""):
			qt.QMessageBox.about(self, "Error!", "Please complete the form and try again.")
			return

		# Catalogue path
		self.parent.catalogue_path = self.catBox.text()
		if not os.path.exists(self.parent.catalogue_path):
			qt.QMessageBox.about(self, "Error!", "You must provide a valid catalogue path.")
			return

		# Catalogue name
		self.parent.catalogue_name = self.catName.text()
		if not type(self.parent.catalogue_name) == str:
			qt.QMessageBox.about(self, "Error!", "The catalogue name must be a String.")
			return

		# Archive path
		self.parent.archive_path = self.arcBox.text()
		if self.parent.archive_path == "":
			pass
		elif not os.path.exists(self.parent.archive_path):
			qt.QMessageBox.about(self, "Error!", "You must provide a valid archive path.")
			return

		# Archive format
		self.parent.archive_format = self.arcFormBox.text()
		
		# Receiver file
		self.parent.receiver_file  = self.recBox.text()
		if not os.path.exists(self.parent.receiver_file):
			qt.QMessageBox.about(self, "Error!", "You must provide a valid receiver file.")
			return

		if self.tele_radio.isChecked():
			self.parent.catalogue_type = "teleseismic"
		elif self.local_radio.isChecked():
			self.parent.catalogue_type = "local"
		else:
			qt.QMessageBox.about(self, "Error!", "You must select a catalogue type.")
			return

		# Send accept signal to Dialog
		self.accept()

	def actionReject(self):
		# Send reject signal to Dialog
		self.reject()

class TelePhaseDialogue(qt.QDialog):

	def __init__(self):
		super(TelePhaseDialogue, self).__init__()

		self.phases = []

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/tele_phase_dialogue.ui', self)

		# Connect to the add phase button
		self.button_addPhase.clicked.connect(self.addPhase)

		# Populate the station list
		self.model = QtGui.QStandardItemModel(self.phase_list)

		# Turn off editing triggers
		self.phase_list.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)

		self.setWindowTitle('PySplit - Teleseismic phase selection')
		self.show()

	def addPhase(self):
		# Grab the current phase to be added
		phase = self.input_phaseType.currentText()

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
			self.phase_list.setModel(self.model)

class LocalInputDialogue(qt.QDialog):

	def __init__(self, parent):
		super(LocalInputDialogue, self).__init__()

		self.parent = parent

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/local_input_dialogue.ui', self)

		# Display the default input file option (SeisLoc)
		self.inputBox.setCurrentIndex(1)
		self.input_type = "SeisLoc"

		# Connect to all actions
		self.inpButton1.clicked.connect(self.browseLocalFile)
		self.inpButton2.clicked.connect(self.browseLocalPath)
		self.buttonBox.accepted.connect(self.actionAccept)
		self.buttonBox.rejected.connect(self.actionReject)
		self.inputType.currentIndexChanged.connect(self.inputSelect)

		self.setWindowTitle('PySplit - Local catalogue input type')
		self.show()

	def browseLocalFile(self):

		filename = qt.QFileDialog.getOpenFileName(self, 'Open file')

		self.input_file = filename[0]

		self.fileInput.setText(self.input_file)

	def browseLocalPath(self):

		pathname = qt.QFileDialog.getExistingDirectory(self, 'Choose SeisLoc directory')

		self.pathInput.setText(pathname)

	def actionAccept(self):
		# If the inputs are accepted, set the parameter within the parent class
		# Perform checks on the input provided
		if self.input_type == "SeisLoc":
			if not os.path.exists(self.pathInput.text()):
				qt.QMessageBox.about(self, "Error!", "You must provide a valid input directory.")
				return

			self.parent.local_input = self.pathInput.text()
			self.parent.data_source = "SeisLoc"
			self.parent.picked = self.picked_radio.isChecked()

		else:
			if not os.path.isfile(self.fileInput.text()):
				qt.QMessageBox.about(self, "Error!", "You must provide a valid input file.")
				return

			# After all checks, set the parameter within the parent class
			self.parent.local_input = self.fileInput.text()
			self.parent.data_source = self.inputType.currentText()
			self.parent.picked = self.picked_Radio.isChecked()

		# Send accept signal to Dialog
		self.accept()

	def actionReject(self):
		# Send reject signal to Dialog
		self.reject()

	def inputSelect(self):
		# Get the current index
		self.input_type = self.inputType.currentText()

		# Handle each different type of input
		if self.input_type == "SeisLoc":
			self.inputBox.setCurrentIndex(1)
		else:
			self.inputBox.setCurrentIndex(0)


class DefaultFilterDialogue(qt.QDialog):

	def __init__(self):
		super(DefaultFilterDialogue, self).__init__()

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/default_filter_dialogue.ui', self)

		self.setWindowTitle('PySplit - Set default filter')
		self.show()

class CustomPickDialogue(qt.QDialog):

	def __init__(self):
		super(CustomPickDialogue, self).__init__()

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/custom_phase_dialogue.ui', self)

		self.setWindowTitle('PySplit - Set custom phase pick')
		self.show()

class WadatiWindow(qt.QMainWindow):

	def __init__(self, parent):
		super(WadatiWindow, self).__init__()

		self.ptravels = []
		self.stravels = []
		self.sptimes  = []
		self.stations = []
		self.dists    = []

		self.parent = parent

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/wadati_window.ui', self)

		# Set plot details (axes labels etc)
		self.wadatiPlot.canvas.ax.set_xlabel("P traveltime / s", fontsize=10)
		self.wadatiPlot.canvas.ax.set_ylabel("S - P traveltime / s", fontsize=10)

		self.connect()

		self.setWindowTitle('PySplit - Wadati plot window')
		self.show()

	def connect(self):
		# Calculate Vp/Vs button
		self.button_calcVpVs.clicked.connect(self.calcVpVs)

		# Plot connection
		self.wadatiPlot.canvas.mpl_connect('pick_event', self._onPick)
		self.distancePlot.canvas.mpl_connect('pick_event', self._onPick)

	def calcVpVs(self):
		# Insert code to calc best fitting line to p travel and s-p times
		# Should be a straight line
		pass

	def addPick(self, ptravel, stravel, dist, station):
		# Calculate sptime
		sptime = stravel - ptravel
		# Add a p traveltime and SP-time to the plot
		self.ptravels.append(ptravel)
		self.stravels.append(stravel)
		self.sptimes.append(sptime)
		self.dists.append(dist)
		self.stations.append(station)

		self.plotWadati(station)
		self.plotDistance(station)

	def _onPick(self, event):
	    artist = event.artist
	    label = artist.get_label()

	    self.label_station.setText(label)

	def plotWadati(self, station):
		wadati_canvas = self.wadatiPlot.canvas

		# Clear the canvas
		wadati_canvas.ax.clear()

		# Set plot details (axes labels etc)
		wadati_canvas.ax.set_xlabel("P traveltime / s", fontsize=10)
		wadati_canvas.ax.set_ylabel("S - P traveltime / s", fontsize=10)

		# Try rescaling the image now
		wadati_canvas.ax.set_aspect('auto')

		tolerance = 10
		for i in range(len(self.ptravels)):
			wadati_canvas.ax.scatter(self.ptravels[i], self.sptimes[i], 12, marker='o', color='k', picker=tolerance, zorder=10, label=self.stations[i])

		wadati_canvas.draw_idle()

	def plotDistance(self, station):
		distance_canvas = self.distancePlot.canvas

		# Clear the canvas
		distance_canvas.ax.clear()

		# Set plot details (axes labels etc)
		distance_canvas.ax.set_xlabel("Epicentral distance / km", fontsize=10)
		distance_canvas.ax.set_ylabel("Traveltime / s", fontsize=10)

		# Rescale the image now
		distance_canvas.ax.set_aspect('auto')

		tolerance = 10
		for i in range(len(self.ptravels)):
			distance_canvas.ax.scatter(self.dists[i], self.ptravels[i], 12, marker='o', color="red", picker=tolerance, zorder=10, label=self.stations[i])
			distance_canvas.ax.scatter(self.dists[i], self.stravels[i], 12, marker='o', color="blue", picker=tolerance, zorder=10, label=self.stations[i])

		distance_canvas.draw_idle()

class PickingWindow(qt.QMainWindow):

	def __init__(self, catalogue, catalogue_name, filt=None, event=None, station=None):
		super(PickingWindow, self).__init__()

		# Parse args
		self.catalogue = catalogue
		self.catalogue_path = self.catalogue.path
		self.catalogue_name = catalogue_name
		self.event   = event
		self.station = station

		# Set default filter
		self.default_filter = filt
		self.filt           = filt

		# Set default pick type and line color
		self.pick_type = "P"
		self.pick_line_color = "red"
		self.pick_polarity = ""

		# Initialise pick tracker
		self.pick_lines = {}

		# Initialise trackers for whether or not multiple events/stations are being picked on
		self.evts  = False
		self.stats = False

		# Initialise trace rejection tracker
		self.trace_removed = False

		# Initialise toggle trackers
		self.ctrl_toggle  = False
		self.shift_toggle = False

		# If picking all stations for a given event
		if station == None and event != None:
			self.stats = True

		# If picking all events at a given station
		elif event == None and station != None:
			self.evts = True

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/trace_window.ui', self)

		# If just picking a single event at a single station
		if not self.evts and not self.stats:
			# Disable the next trace function
			self.button_nextTrace.setEnabled(False)
			self.button_lastTrace.setEnabled(False)

			self._updateEventInformation(self.event)
			print(self.event)
			self._updateStationInformation(self.station)

		# If picking all events at a given station
		elif self.evts:
			# Load the events
			station_path = "{}/data/{}".format(self.catalogue_path, self.station.upper())
			events = glob.glob('{}/*.z'.format(station_path))
			self.events = []
			for event in events:
				head, tail  = os.path.split(event)
				self.events.append(tail.split(".")[1])

			self.counter = 0

			try:
				self.event = self.events[self.counter]
				print(self.event)

				self._updateStationInformation(self.station)
				self._updateEventInformation(self.event)
			except IndexError:
				qt.QMessageBox.about(self, "Error!", "There are no events with recorded arrivals at this station.")
				return

		# If picking all stations for a given event
		elif self.stats:
			# Load stations
			event_path = "{}/data/*/event.{}.*.z".format(self.catalogue_path, self.event)
			stations = glob.glob(event_path)
			self.stations = []
			for station in stations:
				head, tail = os.path.split(station)
				self.stations.append(tail.split(".")[2])

			self.counter = 0

			try:
				self.station = self.stations[self.counter]

				self._updateStationInformation(self.station)
				self._updateEventInformation(self.event)
			except IndexError:
				qt.QMessageBox.about(self, "Error!", "There are no stations with recorded arrivals for this event.")
				return

			# Open up Wadati plot
			self.wadatiWindow = WadatiWindow(self)

		if self.default_filter != None:
			self.input_filtType.setCurrentText(self.default_filter["filt_type"])
			self.input_noPoles.setCurrentText(str(self.default_filter["no_poles"]))
			if self.default_filter["zerophase"]:
				self.input_zerophase.setChecked(True)
			self.input_lowFreq.setText(str(self.default_filter["low_freq"]))
			self.input_highFreq.setText(str(self.default_filter["high_freq"]))

		self.plotEvent(self.event, self.station)

		self.connect()

		self.setWindowTitle('PySplit - Shear Wave Splitting software')
		self.show()

	def connect(self):
		# File menu connections
		self.actionSave.triggered.connect(self.saveTrace)

		# Settings menu connections
		self.actionDefault_filter.triggered.connect(self.defaultFilter)
		self.actionP_pick.triggered.connect(self.pPick)
		self.actionS_pick.triggered.connect(self.sPick)
		self.actionCustom_pick.triggered.connect(self.customPick)

		# Connect to filter buttons
		self.cidafilter = self.button_applyFilter.clicked.connect(self.applyFilter)
		self.cidrfilter = self.button_removeFilter.clicked.connect(self.removeFilter)

		# Connect to plot option buttons
		self.cidnexttrace = self.button_nextTrace.clicked.connect(self.nextTrace)
		self.cidlasttrace = self.button_lastTrace.clicked.connect(self.previousTrace)
		self.cidrjcttrace = self.button_rejectTrace.clicked.connect(self.rejectTrace)
		self.cidlastzoom  = self.button_toggleLims.clicked.connect(self.toggleLims)
		self.cidresetplot = self.button_resetPlot.clicked.connect(self.resetPlot)

	def saveTrace(self):
		if self.evt != None:
			print("Saving event...")
			self.evt.save_event("{}/picks".format(self.catalogue_path))
		else:
			return

	def rejectTrace(self):
		if self.evt != None:
			print("Deleting trace...")
			traces = glob.glob(self.evt.file_path)
			for trace in traces:
				os.remove(trace)

			# Pass information back to parent and remove arrival from available arrivals?

			# Delete event or station from events/stations list
			if self.evts:
				del self.events[self.counter]

			if self.stats:
				del self.stations[self.counter]

			self.trace_removed = True

			self.nextTrace()
		else:
			return

	def plotconnect(self):
		# Connect each canvas to track motion - this will create a 
		self.cidzcompmove = self.z_comp_plot.canvas.mpl_connect('motion_notify_event', self._onMove)
		self.cidncompmove = self.n_comp_plot.canvas.mpl_connect('motion_notify_event', self._onMove)
		self.cidecompmove = self.e_comp_plot.canvas.mpl_connect('motion_notify_event', self._onMove)

		# Connect to clicks
		self.cidzcompclick = self.z_comp_plot.canvas.mpl_connect('button_press_event', self._onClick)
		self.cidncompclick = self.n_comp_plot.canvas.mpl_connect('button_press_event', self._onClick)
		self.cidecompclick = self.e_comp_plot.canvas.mpl_connect('button_press_event', self._onClick)

		# Connect to releases
		self.cidzcomprelease = self.z_comp_plot.canvas.mpl_connect('button_release_event', self._onRelease)
		self.cidncomprelease = self.n_comp_plot.canvas.mpl_connect('button_release_event', self._onRelease)
		self.cidecomprelease = self.e_comp_plot.canvas.mpl_connect('button_release_event', self._onRelease)

	def plotdisconnect(self):
		# Disconnect each canvas to track motion - this will create a 
		self.z_comp_plot.canvas.mpl_disconnect(self.cidzcompmove)
		self.n_comp_plot.canvas.mpl_disconnect(self.cidncompmove)
		self.e_comp_plot.canvas.mpl_disconnect(self.cidecompmove)

		# Disconnect to clicks
		self.z_comp_plot.canvas.mpl_disconnect(self.cidzcompclick)
		self.n_comp_plot.canvas.mpl_disconnect(self.cidncompclick)
		self.e_comp_plot.canvas.mpl_disconnect(self.cidecompclick)

		# Disconnect to releases
		self.z_comp_plot.canvas.mpl_disconnect(self.cidzcomprelease)
		self.n_comp_plot.canvas.mpl_disconnect(self.cidncomprelease)
		self.e_comp_plot.canvas.mpl_disconnect(self.cidecomprelease)

	def keyPressEvent(self, event):
		# Toggle Ctrl Modifier on
		if event.key() == Qt.Key_Control:
			self.ctrl_toggle = True

		# Toggle Shift Modifier on
		if event.key() == Qt.Key_Shift:
			self.shift_toggle = True

		# Up polarity indicator
		if event.key() == Qt.Key_U:
			print("U")
			if not self.pick_line:
				return
			else:
				# Set pick polarity
				self.pick_polarity = "U"

				# Update label
				self.label_pickPolarity.setText(self.pick_polarity)

				# Add to event
				self.evt._add_stat(stat="polarity", value=self.pick_polarity, pick_type=self.pick_type)

		# Down polarity indicator
		if event.key() == Qt.Key_D:
			print("D")
			if not self.pick_line:
				return
			else:
				# Set pick polarity
				self.pick_polarity = "D"

				# Update label
				self.label_pickPolarity.setText(self.pick_polarity)

	def keyReleaseEvent(self, event):
		# Toggle Ctrl Modifier off
		if event.key() == Qt.Key_Control:
			self.ctrl_toggle = False

		# Toggle Shift Modifier off
		if event.key() == Qt.Key_Shift:
			self.shift_toggle = False

	def pPick(self):
		self.pick_type = "P"
		self.pick_line_color = "red"

		# Reset polarity text
		self.label_pickPolarity.setText("")

	def sPick(self):
		self.pick_type = "S"
		self.pick_line_color = "blue"

		# Reset polarity text
		self.label_pickPolarity.setText("")

	def customPick(self):
		# Only run if the radio button is being toggled on.
		if not self.c_radio.isChecked():
			return

		# Reset polarity text
		self.label_pickPolarity.setText("")

		# Open custom phase dialogue
		self.customPickDialogue = CustomPickDialogue()

		if self.customPickDialogue.exec_():
			# Read the custom pick type
			try:
				self.pick_type = self.customPickDialogue.input_pickType.currentText()
				self.label_customPhase.setText(self.pick_type)

				self.pick_line_color = "orange"
			except ValueError:
				qt.QMessageBox.about(self, "Error!", "You need to specify a phase to pick!")
		else:
			return

	def defaultFilter(self):
		self.defaultFilterDialogue = DefaultFilterDialogue()

		if self.defaultFilterDialogue.exec_():
			# Read the default filter parameters
			try:
				filt_type = self.defaultFilterDialogue.input_filtType.currentText()
				no_poles  = int(self.defaultFilterDialogue.input_noPoles.currentText())
				zerophase = self.defaultFilterDialogue.input_zerophase.isChecked()
				low_freq  = float(self.defaultFilterDialogue.input_lowFreq.text())
				high_freq = float(self.defaultFilterDialogue.input_highFreq.text())

				self.filt = {'filt_type': filt_type,
							 'no_poles': no_poles,
							 'zerophase': zerophase,
							 'low_freq': low_freq,
							 'high_freq': high_freq
							}
			except ValueError:
				qt.QMessageBox.about(self, "Error!", "You appear to have tried to use an incomplete/incorrect filter. Fill in all of the options and try again.")
		else:
			return


	def _onMove(self, event):
		# Temp variables for accessing canvases
		z_canvas = self.z_comp_plot.canvas
		n_canvas = self.n_comp_plot.canvas
		e_canvas = self.e_comp_plot.canvas

		if event.inaxes is z_canvas.ax or n_canvas.ax or e_canvas.ax:
			if self.shift_toggle:
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
				self.zcursor.set_data([x, x], self.zydat)
				self.ncursor.set_data([x, x], self.nydat)
				self.ecursor.set_data([x, x], self.eydat)

			# restore the background region
			z_canvas.restore_region(self.z_background)
			n_canvas.restore_region(self.n_background)
			e_canvas.restore_region(self.e_background)

			if self.shift_toggle:
				# Set rectangle values
				z_trace_rectangle = Rectangle((xpress, ypress), dx, dy,
										  edgecolor='red', fill=False)
				n_trace_rectangle = Rectangle((xpress, ypress), dx, dy,
										  edgecolor='red', fill=False)
				e_trace_rectangle = Rectangle((xpress, ypress), dx, dy,
										  edgecolor='red', fill=False)
				# Add the rectangle patch to each canvas
				z_canvas.ax.add_patch(z_trace_rectangle)
				n_canvas.ax.add_patch(n_trace_rectangle)
				e_canvas.ax.add_patch(e_trace_rectangle)

				# Redraw just the rectangle
				z_canvas.ax.draw_artist(z_trace_rectangle)
				n_canvas.ax.draw_artist(n_trace_rectangle)
				e_canvas.ax.draw_artist(e_trace_rectangle)

			else:
				# Redraw just the current cursor
				z_canvas.ax.draw_artist(self.zcursor)
				n_canvas.ax.draw_artist(self.ncursor)
				e_canvas.ax.draw_artist(self.ecursor)

			# Plot any lines that are currently being stored
			self._replotLines()				

			# blit just the redrawn area
			z_canvas.blit(z_canvas.ax.bbox)
			n_canvas.blit(n_canvas.ax.bbox)
			e_canvas.blit(e_canvas.ax.bbox)

	def _onClick(self, event):
		# Check for shift modifier - this enables the draggable zoom
		if self.shift_toggle:
			if self._trace_drag_lock is not None:
				return

			self._zoom_click = True

			xpress, ypress = event.xdata, event.ydata
			self._trace_drag_lock = self

			# Grab the boundary values
			self.trace_click = xpress, ypress

		else:
			z_canvas = self.z_comp_plot.canvas
			n_canvas = self.n_comp_plot.canvas
			e_canvas = self.e_comp_plot.canvas

			dt = self.evt.Z_comp.stats.delta

			adjusted_xdata = round(event.xdata / dt) * dt

			# Left-clicking handles the window start time
			if event.button == 1 and not self.ctrl_toggle:
				# Set the window start line to be redrawn on move
				self.w_beg_line = True

				# Add the window beginning to the event stats
				self.evt._add_stat(stat="window_beg", value=adjusted_xdata)

				# Make a vertical line artist
				self.z_window_beg = z_canvas.ax.axvline(adjusted_xdata, linewidth=1, color="green", animated=True)
				self.n_window_beg = n_canvas.ax.axvline(adjusted_xdata, linewidth=1, color="green", animated=True)
				self.e_window_beg = e_canvas.ax.axvline(adjusted_xdata, linewidth=1, color="green", animated=True)

				# Restore the background region
				z_canvas.restore_region(self.z_background)
				n_canvas.restore_region(self.n_background)
				e_canvas.restore_region(self.e_background)

				# Plot any lines that are currently being stored
				self._replotLines()

				# blit the redrawn area
				z_canvas.blit(z_canvas.ax.bbox)
				n_canvas.blit(n_canvas.ax.bbox)
				e_canvas.blit(e_canvas.ax.bbox)

			# Middle-clicking handles the arrival pick time
			if event.button == 2 or (event.button == 1 and self.ctrl_toggle):
				# Set the pick line to be redrawn on move
				self.pick_line = True

				# Make a UTCDateTime for the pick time
				pick_time = self.evt.starttime + adjusted_xdata

				# Add the pick to the event stats
				self.evt._add_stat(stat="pick", value=adjusted_xdata, pick_type=self.pick_type)

				# Set pick time label
				self.label_pickTime.setText(pick_time.isoformat())
				self.label_pickPhase.setText(self.pick_type)

				# Make a vertical line artist
				self.z_pick = z_canvas.ax.axvline(adjusted_xdata, linewidth=1, color=self.pick_line_color, animated=True)
				self.n_pick = n_canvas.ax.axvline(adjusted_xdata, linewidth=1, color=self.pick_line_color, animated=True)
				self.e_pick = e_canvas.ax.axvline(adjusted_xdata, linewidth=1, color=self.pick_line_color, animated=True)

				self.pick_lines[self.pick_type] = [self.z_pick, self.n_pick, self.e_pick]

				# Restore the background region
				z_canvas.restore_region(self.z_background)
				n_canvas.restore_region(self.n_background)
				e_canvas.restore_region(self.e_background)

				# Plot any lines that are currently being stored
				self._replotLines()

				# blit the redrawn area
				z_canvas.blit(z_canvas.ax.bbox)
				n_canvas.blit(n_canvas.ax.bbox)
				e_canvas.blit(e_canvas.ax.bbox)

			# Right-clicking handles the window end time
			if event.button == 3:
				# Set the window end line to be redrawn on move
				self.w_end_line = True

				# Add the window ending to the event stats
				self.evt._add_stat(stat="window_end", value=adjusted_xdata)

				# Make a vertical line artist
				self.z_window_end = z_canvas.ax.axvline(adjusted_xdata, linewidth=1, color="green", animated=True)
				self.n_window_end = n_canvas.ax.axvline(adjusted_xdata, linewidth=1, color="green", animated=True)
				self.e_window_end = e_canvas.ax.axvline(adjusted_xdata, linewidth=1, color="green", animated=True)

				# Restore the background region
				z_canvas.restore_region(self.z_background)
				n_canvas.restore_region(self.n_background)
				e_canvas.restore_region(self.e_background)

				# Plot any lines that are currently being stored
				self._replotLines()

				# blit the redrawn area
				z_canvas.blit(z_canvas.ax.bbox)
				n_canvas.blit(n_canvas.ax.bbox)
				e_canvas.blit(e_canvas.ax.bbox)

	def _replotLines(self):
		# For clarity
		z_canvas = self.z_comp_plot.canvas
		n_canvas = self.n_comp_plot.canvas
		e_canvas = self.e_comp_plot.canvas

		if self.w_beg_line:
			# Redraw the window line
			z_canvas.ax.draw_artist(self.z_window_beg)
			n_canvas.ax.draw_artist(self.n_window_beg)
			e_canvas.ax.draw_artist(self.e_window_beg)

		if self.pick_line:
			# Redraw the pick lines
			for phase, pick_lines in self.pick_lines.items():
				z_canvas.ax.draw_artist(pick_lines[0])
				n_canvas.ax.draw_artist(pick_lines[1])
				e_canvas.ax.draw_artist(pick_lines[2])

		if self.w_end_line:
			# Redraw the window line
			z_canvas.ax.draw_artist(self.z_window_end)
			n_canvas.ax.draw_artist(self.n_window_end)
			e_canvas.ax.draw_artist(self.e_window_end)

	def _onRelease(self, event):
		if self.shift_toggle:
		
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

			self.plotEvent(self.event, self.station, replot=True)

	def applyFilter(self):
		# Read the default filter parameters
		try:
			filt_type = self.input_filtType.currentText()
			no_poles  = int(self.input_noPoles.currentText())
			zerophase = self.input_zerophase.isChecked()
			low_freq  = float(self.input_lowFreq.text())
			high_freq = float(self.input_highFreq.text())

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
		self.plotEvent(self.event, self.station, replot=True)

	def removeFilter(self):
		self.plotdisconnect()
		self.filt = None
		self.evt.remove_filter()
		self.plotEvent(self.event, self.station, replot=True)

	def toggleLims(self):
		# Overwrite limits
		self.lims, self.previous_lims = self.previous_lims, self.lims

		# Replot the traces
		self.plotEvent(self.event, self.station, replot=True)

	def nextTrace(self):
		# If previous trace was removed, do not advance counter
		if not self.trace_removed:
			# Advance counter
			self.counter += 1

		# Check if P and S have been picked (and window exists)
		if "P_manual" in self.evt.picks.keys() and "S_manual" in self.evt.picks.keys() and self.stats:
			ptravel = self.evt.starttime - self.evt.otime
			ptravel += self.evt.picks["P_manual"]["rtime"]
			stravel = self.evt.starttime - self.evt.otime
			stravel += self.evt.picks["S_manual"]["rtime"]
			self.wadatiWindow.addPick(ptravel, stravel, self.evt.distance, self.station)

		# Save any picks
		self.saveTrace()

		# Reset to P pick
		self.p_radio.setChecked(True)

		# Reset trace removal tracker
		self.trace_removed = False

		# Reset pick trackers
		self.pick_lines = {}
		self.pick_polarity = ""
		self.pick_time = ""
		self.pick_phase = ""

		# Set all the labels
		self.label_pickPolarity.setText(self.pick_polarity)
		self.label_pickTime.setText(self.pick_time)
		self.label_pickPhase.setText(self.pick_phase)

		if self.evts:
			# Check not on last event
			if self.counter == len(self.events):
				qt.QMessageBox.about(self, "Warning!", "Looped round to first event.")
				self.counter = 0

			# Change current event
			self.event = self.events[self.counter]

			# Update event information
			self._updateEventInformation(self.event)

		if self.stats:
			# Check not on last station
			if self.counter == len(self.stations):
				qt.QMessageBox.about(self, "Warning!", "Looped round to first station.")
				self.counter = 0

			# Change current station
			self.station = self.stations[self.counter]

			# Update station information
			self._updateStationInformation(self.station)

		self.plotdisconnect()
		self.plotEvent(self.event, self.station)

	def previousTrace(self):
		# Reverse counter
		self.counter += -1

		# Check if P and S have been picked (and window exists)
		if "P_manual" in self.evt.picks.keys() and "S_manual" in self.evt.picks.keys() and self.stats:
			ptravel = self.evt.starttime - self.evt.otime
			ptravel += self.evt.picks["P_manual"]["rtime"]
			stravel = self.evt.starttime - self.evt.otime
			stravel += self.evt.picks["S_manual"]["rtime"]
			self.wadatiWindow.addPick(ptravel, stravel, self.station)

		# Save any picks
		self.saveTrace()

		# Reset to P pick
		self.p_radio.setChecked(True)

		# Reset trace removal tracker
		self.trace_removed = False

		# Reset pick tracker
		self.pick_lines = {}
		self.pick_polarity = ""
		self.pick_time = ""
		self.pick_phase = ""

		# Set all the labels
		self.label_pickPolarity.setText(self.pick_polarity)
		self.label_pickTime.setText(self.pick_time)
		self.label_pickPhase.setText(self.pick_phase)

		if self.evts:
			# Check not trying to go below 0th position
			if self.counter == -1:
				qt.QMessageBox.about(self, "Warning!", "Looped round to last event.")
				self.counter = len(events) - 1

			# Change current event
			self.event = self.events[self.counter]

			# Update event information
			self._updateEventInformation(self.event)

		if self.stats:
			# Check not trying to go below 0th position
			if self.counter == -1:
				qt.QMessageBox.about(self, "Warning!", "Looped round to last station.")
				self.counter = len(stations) - 1

			# Change current station
			self.station = self.stations[self.counter]

			# Update station information
			self._updateStationInformation(self.station)

		self.plotdisconnect()
		self.plotEvent(self.event, self.station)

	def resetPlot(self):
		self.plotdisconnect()
		self.lims = None
		self.filt = self.default_filter
		self.plotEvent(self.event, self.station)

	def plotEvent(self, event, station, replot=False):

		z_canvas = self.z_comp_plot.canvas
		n_canvas = self.n_comp_plot.canvas
		e_canvas = self.e_comp_plot.canvas

		# Clear the canvases
		z_canvas.ax.clear()
		n_canvas.ax.clear()
		e_canvas.ax.clear()

		# Create background variables
		self.z_background = None
		self.n_background = None
		self.e_background = None

		# Create the lock and click variables
		self._zoom_click      = False
		self._trace_drag_lock = None
		self.trace_click      = None

		if not replot:
			# Add limits variable
			self.previous_lims = None
			self.lims = None

			# Set filter to default filter
			self.filt = self.default_filter

			# Make the window and pick lines not show on move
			self.w_beg_line = False
			self.pick_line  = False
			self.w_end_line = False

			# Create an instance of the Event class
			self.evt = evt.Event("{}/data/{}/event.{}.{}.*".format(self.catalogue_path, station.upper(), event, station.upper()), self.event_info, self.station_info)

			# Look up any picks
			# try:
			# 	pick_files = glob.glob("{}/picks/{}/event.{}.{}*".format(self.catalogue_path, station.upper(), event, station.upper()))
			# 	for pick_file in pick_files:

			# 		# Handle window start files
			# 		if pick_file[-2:] == "wb":
			# 			# Read in window start time
			# 			window = pd.read_csv(pick_file, sep=" ", header=None)
			# 			w_beg_time = window.iloc[0][0]

			# 			# Make a vertical line artist
			# 			z_canvas.ax.axvline(w_beg_time, linewidth=1, color="green")
			# 			n_canvas.ax.axvline(w_beg_time, linewidth=1, color="green")
			# 			e_canvas.ax.axvline(w_beg_time, linewidth=1, color="green")

			# 		# Handle window end files
			# 		elif pick_file[-2:] == "we":
			# 			# Read in window end time
			# 			window = pd.read_csv(pick_file, sep=" ", header=None)
			# 			w_end_time = window.iloc[0][0]

			# 			# Make a vertical line artist
			# 			z_canvas.ax.axvline(w_end_time, linewidth=1, color="green")
			# 			n_canvas.ax.axvline(w_end_time, linewidth=1, color="green")
			# 			e_canvas.ax.axvline(w_end_time, linewidth=1, color="green")

			# 		# Handle picks
			# 		else:
			# 			pick_phase = (pick_file.split(".")[-1]).split("_")[0]

			# 			pick = pd.read_csv(pick_file, sep=" ", header=None)

			# 			pick_time = UTCDateTime("{}T{}:{}".format(pick.iloc[0][0], pick.iloc[0][1], pick.iloc[0][2])) - self.evt.starttime

			# 			# Handle P picks
			# 			if pick_phase == "P":
			# 				self.pick_line_color = "red"

			# 			# Handle S picks
			# 			elif pick_phase == "S":
			# 				self.pick_line_color = "blue"

			# 			# Handle other picks
			# 			else:
			# 				self.pick_line_color = "orange"

			# 			z_canvas.ax.axvline(pick_time, linewidth=1, color=self.pick_line_color)
			# 			n_canvas.ax.axvline(pick_time, linewidth=1, color=self.pick_line_color)
			# 			e_canvas.ax.axvline(pick_time, linewidth=1, color=self.pick_line_color)

			# 	self.pick_line_color = "red"

			# except FileNotFoundError:
			# 	print("No picks exist for this station")
			# 	pass

		else:
			self._replotLines()

		# Check if a filter has been specified, and apply it if so
		if not self.filt == None:
			self.evt.filter_obspy(filt_type=self.filt["filt_type"], minfreq=self.filt["low_freq"], maxfreq=self.filt["high_freq"], n_poles=self.filt["no_poles"], zero_phase=self.filt["zerophase"])

		# Plot the traces
		self.evt.plot_traces(z_canvas.ax, n_canvas.ax, e_canvas.ax, lims=self.lims)

		# Connect to each trace to grab the background once Qt has done resizing
		z_canvas.mpl_connect('draw_event', self._zDrawEvent)
		n_canvas.mpl_connect('draw_event', self._nDrawEvent)
		e_canvas.mpl_connect('draw_event', self._eDrawEvent)

		# Draw when Qt has done all resizing
		z_canvas.draw_idle()
		n_canvas.draw_idle()
		e_canvas.draw_idle()

		# Initialise the cursor to track mouse position on the axes
		self.zcursor = z_canvas.ax.axvline(5, linewidth=1, color='0.5', animated=True)
		_, self.zydat = self.zcursor.get_data()
		self.ncursor = n_canvas.ax.axvline(5, linewidth=1, color='0.5', animated=True)
		_, self.nydat = self.ncursor.get_data()
		self.ecursor = e_canvas.ax.axvline(5, linewidth=1, color='0.5', animated=True)
		_, self.eydat = self.ecursor.get_data()

		self.plotconnect()

	def _zDrawEvent(self, evt):
		self.z_background = evt.canvas.copy_from_bbox(evt.canvas.ax.bbox)

	def _nDrawEvent(self, evt):
		self.n_background = evt.canvas.copy_from_bbox(evt.canvas.ax.bbox)

	def _eDrawEvent(self, evt):
		self.e_background = evt.canvas.copy_from_bbox(evt.canvas.ax.bbox)

	def _updateStationInformation(self, station):
		station_info = self.catalogue.receiver_df.query('stat == @station')
		self.station_info = station_info
		self.label_statName.setText(station)
		self.label_statLon.setText(f"{station_info.lon.values[0]:.4f}")
		self.label_statLat.setText(f"{station_info.lat.values[0]:.4f}")
		self.label_statElev.setText(f"{station_info.dep.values[0]:.4f}")
		self.label_statDep.setText(station_info.st_dep.values[0].isoformat().split("T")[0])
		self.label_statRet.setText(station_info.et_dep.values[0].isoformat().split("T")[0])

	def _updateEventInformation(self, event):
		event_info = self.catalogue.source_df.query('sourceid == @event')
		self.event_info = event_info
		self.label_eventOdate.setText(event_info.otime.values[0].isoformat().split("T")[0])
		self.label_eventOtime.setText(event_info.otime.values[0].isoformat().split("T")[1])
		self.label_eventLon.setText(f"{event_info.evlon.values[0]:.4f}")
		self.label_eventLat.setText(f"{event_info.evlat.values[0]:.4f}")
		self.label_eventDepth.setText(f"{event_info.evdep.values[0]:.4f}")
		if type(event_info.evmag.values[0]) == float:
			self.label_eventMag.setText(f"{event_info.evmag.values[0]:.2f}")
		else:
			self.label_eventMag.setText(event_info.evmag.values[0])
		self.label_eventID.setText(str(event_info.sourceid.values[0]))


if __name__ == "__main__":

	app = qt.QApplication(sys.argv)
	pysplit = PySplit()
	sys.exit(app.exec_())