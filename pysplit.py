#!/usr/bin/env python3

"""
The Map class - generate an instance of Map and plot the geographic
data. The data points are interactive, and clicking one will create
an instance of the Picker class, which will allow the window to be
manually selected.

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
		uic.loadUi('ui_files/pysplit.ui', self)

		# Connect to actions and interactive widgets
		self.connect()

		# Hide the map widget initially
		self.mpl.hide()
		self.teleseismic_input_page.hide()
		self.status_page.hide()
		self.catalogue_info_page.hide()
		self.station_info_page.hide()

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

		# Station list connection
		self.station_list.doubleClicked.connect(self.stationSelect)

		# Event list connection
		self.event_list.doubleClicked.connect(self.eventSelect)

		# Catalogue generation button connection
		self.tele_load_button.clicked.connect(self.generate_teleseismic_catalogue)

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

	def keyPressEvent(self, event):
		if event.key() == Qt.Key_Control:
			print("Ctrl toggled on")
			self.ctrl_toggle = True
		if event.key() == Qt.Key_Shift:
			print("Shift toggled on")
			self.shift_toggle = True

	def keyReleaseEvent(self, event):
		if event.key() == Qt.Key_Control:
			print("Ctrl toggled off")
			self.ctrl_toggle = False
		if event.key() == Qt.Key_Shift:
			print("Shift toggled off")
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
			self.catalogue.station_points.set_visible(True)
		else:
			self.catalogue.station_points.set_visible(False)

	def pickStateChanged(self, int):
		if self.picked_plot_tickbox.isChecked():
			self.catalogue.pick_points.set_visible(True)
		else:
			self.catalogue.pick_points.set_visible(False)

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

	    # If event
	    if "EVENT" in label:
	    	evt = label.split(": ")[1]
	    	self._updateEventInformation(evt)

	    	# Open an instance of Picker (TESTING)
	    	if self.actionPick_mode.isChecked():
		    	self.pickerWindow = PickingWindow(self.catalogue)
		    	self.pickerWindow.plot_event(1261, "ASK")

	def defaultFilter(self):
		self.defaultFilterDialogue = DefaultFilterDialogue()

		if self.defaultFilterDialogue.exec_():
			# Read the default filter parameters
			try:
				self.low_freq  = float(self.defaultFilterDialogue.input_lowFreq.text())
				self.high_freq = float(self.defaultFilterDialogue.input_highFreq.text())

				self.filt = [self.low_freq, self.high_freq]
			except:
				qt.QMessageBox.about(self, "Error!", "You did something wrong dickhead, try again")
		else:
			return

	def newCatalogue(self):
		self.newCatalogueDialogue = NewCatalogueDialogue()

		if self.newCatalogueDialogue.exec_():
			# Read in parameters from input dialogue
			self.catalogue_path = self.newCatalogueDialogue.catBox.text()
			self.catalogue_name = self.newCatalogueDialogue.catName.text()
			self.archive_path   = self.newCatalogueDialogue.arcBox.text()
			self.archive_format = self.newCatalogueDialogue.arcFormBox.text()
			self.receiver_file  = self.newCatalogueDialogue.recBox.text()

			if self.newCatalogueDialogue.tele_radio.isChecked():
				self.catalogue_type = "teleseismic"
			elif self.newCatalogueDialogue.local_radio.isChecked():
				self.catalogue_type = "local"
			else:
				qt.QMessageBox.about(self, "Error!", "You must select a catalogue type.")

		else:
			return

		if self.catalogue_type == "teleseismic":
			# Set catalogue type and data source labels
			self.data_source = "IRIS"
			self.label_catTypeDisp.setText("Teleseismic")

			self.status_page.hide()
			self.teleseismic_input_page.show()
			self.station_info_page.hide()
			self.catalogue_info_page.show()

			# Create an instance of the Teleseismic Catalogue class
			self.catalogue = cat.TeleseismicCatalogue("{}/{}".format(self.catalogue_path, self.catalogue_name), self.archive_path, self.receiver_file)

		elif self.catalogue_type == "local":
			# Set catalogue type label
			self.label_catTypeDisp.setText("Local")

			self.teleseismic_input_page.hide()
			self.status_page.show()
			self.catalogue_info_page.hide()
			self.station_info_page.show()

			# Open input file dialogue
			self.getLocalInput()

			# Create an instance of the Local Catalogue class
			self.catalogue = cat.LocalCatalogue("{}/{}".format(self.catalogue_path, self.catalogue_name), self.archive_path, self.receiver_file)

			# Create the catalogue metafile - this contains general information about the catalogue
			# and is read in when an existing catalogue is loaded
			self._generate_catalogue_metafile()
			self.generate_catalogue()

			# Plot the map of the catalogue
			self.plotCatalogueMap()

		# Set the catalogue name label
		self.label_catNameDisp.setText(self.catalogue_name)

		# Populate the station list
		model = QtGui.QStandardItemModel(self.station_list)
		for receiver in self.catalogue.receiver_df.stat.values:
			item = QtGui.QStandardItem(receiver)
			model.appendRow(item)
		self.station_list.setModel(model)
		self.station_list.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)

		# Enable the station and event selection lists
		self.station_list.setEnabled(True)
		self.event_list.setEnabled(True)

	def getLocalInput(self):
		self.localInputDialogue = LocalInputDialogue()

		if self.localInputDialogue.exec_():
			# Read in input file
			self.local_input_file = self.localInputDialogue.inpBox.text()

			if not os.path.isfile(self.local_input_file):
				qt.QMessageBox(self, "Error!", "The input file you have chosen does not exist...")
				return

			# Check if the events have been picked
			if self.localInputDialogue.picked_Radio.isChecked():
				self.picked = True
			else:
				self.picked = False

		else:
			return

		# Set data source
		self.data_source = self.local_input_file


	def plotCatalogueMap(self, replot=False):
		# Transmit a message to the status bar
		self.statusbar.showMessage("Plotting catalogue map...")

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
			self.input_minLon.setText(str(self.catalogue.lon0))
			self.input_maxLon.setText(str(self.catalogue.lon1))
			self.input_minLat.setText(str(self.catalogue.lat0))
			self.input_maxLat.setText(str(self.catalogue.lat1))

		# Draw the canvas
		self.mpl.canvas.draw()

		# Transmit a message to the status bar
		self.statusbar.showMessage("Catalogue load complete.")

	def stationSelect(self, index):
		# Parse station selected
		self.station = self.station_list.model().data(index)

		# Load station information and print to display
		self._updateStationInformation(self.station)

		# Filter and load the corresponding events
		station_path = "{}/{}/data/{}".format(self.catalogue_path, self.catalogue_name, self.station.upper())

		events = glob.glob('{}/*.z'.format(station_path))

		# Populate the event list
		model = QtGui.QStandardItemModel(self.event_list)
		for event in events:
			head, tail = os.path.split(event)
			item = QtGui.QStandardItem("Event {}".format(tail.split(".")[1]))
			model.appendRow(item)
		self.event_list.setModel(model)
		self.event_list.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)

	def eventSelect(self, index):
		# Parse event selected
		self.event = self.event_list.model().data(index).split(" ")[1]

		# Load event information and print to display
		self._updateEventInformation(self.event)

		# If in Pick mode, open up picking window
		if self.actionPick_mode.isChecked():
			self.pickingWindow = PickingWindow(self.catalogue)
			self.pickingWindow.plot_event(self.event, self.station, self.filt)

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
		self.label_eventMag.setText(event_info.evmag.values[0])
		self.label_eventID.setText(str(event_info.sourceid.values[0]))


	def loadCatalogue(self):
		# Open file dialogue and choose catalogue directory
		self.catalogue_path = qt.QFileDialog.getExistingDirectory(self, 'Choose catalogue directory')

		# Parse the catalogue metafile
		self._parse_catalogue_metafile()

		# Set the catalogue name label
		self.label_catNameDisp.setText(self.catalogue_name)
		self.label_catTypeDisp.setText(self.catalogue_type)

		# Create an instance of Catalogue and load in the requisite information
		if self.catalogue_type == "local":
			# Create an instance of Local Catalogue
			self.statusbar.showMessage("Creating an instance of LocalCatalogue...")
			self.catalogue = cat.LocalCatalogue("{}/{}".format(self.catalogue_path, self.catalogue_name), self.archive_path, self.receiver_file)

			# Show the status page catalogues
			self.teleseismic_input_page.hide()
			self.status_page.show()
			self.catalogue_info_page.hide()
			self.station_info_page.show()

		if self.catalogue_type == "teleseismic":
			# Create an instance of Teleseismic Catalogue
			self.statusbar.showMessage("Creating an instance of TeleseismicCatalogue...")
			self.catalogue = cat.TeleseismicCatalogue("{}/{}".format(self.catalogue_path, self.catalogue_name), self.archive_path, self.receiver_file)

			# Show the catalogue info page
			self.status_page.hide()
			self.teleseismic_input_page.show()
			self.station_info_page.hide()
			self.catalogue_info_page.show()

		# Load the catalogue
		self.statusbar.showMessage("Loading catalogue...")
		self.load_events()
		self.load_stations()

		# Populate the station list
		model = QtGui.QStandardItemModel(self.station_list)
		for receiver in self.catalogue.receiver_df.stat.values:
			item = QtGui.QStandardItem(receiver)
			model.appendRow(item)
		self.station_list.setModel(model)
		self.station_list.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)

		# Enable the station and event selection lists
		self.station_list.setEnabled(True)
		self.event_list.setEnabled(True)

		# Plot the map of the catalogue
		self.plotCatalogueMap()

	def generate_teleseismic_catalogue(self):
		pass

	def load_events(self):
		sources = self.catalogue.load_sources()

		if sources:
			self.label_eventStat.setText("Events generated.")
			self.label_eventStat.setStyleSheet("color: rgb(0, 255, 0)")
			self.load_events_button.setEnabled(False)
		else:
			qt.QMessageBox.about(self, "Error!", "Error loading event file.")

	def load_stations(self):
		stations = self.catalogue.load_receivers()

		if stations:
			self.label_statStat.setText("Stations generated.")
			self.label_statStat.setStyleSheet("color: rgb(0, 255, 0)")
			self.load_stats_button.setEnabled(False)
		else:
			qt.QMessageBox.about(self, "Error!", "Error loading station file.")

	def load_arrivals(self):
		arrivals = self.catalogue.load_arrivals()

		if arrivals:
			self.label_arrivalStat.setText("Arrivals generated.")
			self.label_arrivalStat.setStyleSheet("color: rgb(0, 255, 0)")
			self.load_arrivals_button.setEnabled(False)
			self.load_waveforms_button.setEnabled(True)
		else:
			qt.QMessageBox.about(self, "Error!", "The arrivals file does not exist - attempting to generate one...")
			self.catalogue.get_arrivals(self.local_input_file)
			self.load_arrivals()

	def load_waveforms(self):
		try:
			self.catalogue.get_waveforms()
		except:
			qt.QMessageBox.about(self, "Error!", "Some error in there somewhere yo")

		self.label_waveformStat.setText("Waveforms generated.")
		self.label_waveformStat.setStyleSheet("color: rgb(0, 255, 0)")
		self.load_waveforms_button.setEnabled(False)

	def generate_catalogue(self):

		# Get catalogue type
		if self.catalogue_type == "local":
			self.catalogue.generate_catalogue(self.local_input_file)
			#self.catalogue.get_arrivals(self.local_input_file)
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
			self._generate_catalogue_metafile()

			self.teleseismic_input_page.hide()
			self.status_page.show()

			# Plot the map of the catalogue
			self.plotCatalogueMap()

	def _generate_catalogue_metafile(self):
		# Read in the template file
		filein = open('metafile_template.txt')
		src = Template(filein.read())

		# Set the general metadata parameters
		d_cat = {'name': self.catalogue_name,
			 'type': self.catalogue_type,
			 'cpath': self.catalogue_path,
			 'dsource': self.data_source,
			 'cdate': QDateTime.currentDateTime().toUTC().toString(Qt.ISODate),
			 'apath': self.archive_path,
			 'aformat': self.archive_format,
			 'rfile': self.receiver_file,
			 'sdate': self.startDate_input.date().toString(Qt.ISODate),
			 'edate': self.endDate_input.date().toString(Qt.ISODate)}

		# Set the local catalogue parameters
		if self.catalogue_type == "local":
			d_loc = {'lfile': self.local_input_file}

			# Merge the two catalogues
			d = {**d_cat, **d_loc}

		# Set the teleseismic catalogue parameters
		if self.catalogue_type == "teleseismic":
			d_tel = {'minmag': self.minMag_input.text(),
					 'clon': self.centreLon_input.text(),
		 			 'clat': self.centreLat_input.text(),
					 'minrad': self.minRad_input.text(),
		 			 'maxrad': self.maxRad_input.text()}

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
			self.local_input_file = params[10]

		# Parse teleseismic parameters
		if self.catalogue_type == "teleseismic":
			self.minmag = params[11]
			self.clon   = params[12]
			self.clat   = params[13]
			self.minrad = params[14]
			self.maxrad = params[15]

class NewCatalogueDialogue(qt.QDialog):

	def __init__(self):
		super(NewCatalogueDialogue, self).__init__()

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/new_catalogue_dialogue.ui', self)

		self.catDirButton.clicked.connect(self.browseCatalogue)
		self.arcDirButton.clicked.connect(self.browseArchive)
		self.receiverFileButton.clicked.connect(self.browseReceivers)

		self.setWindowTitle('PySplit - New catalogue')
		self.show()

	def browseCatalogue(self):

		self.catalogue_path = qt.QFileDialog.getExistingDirectory(self, 'Choose catalogue directory')

		self.catBox.setText(self.catalogue_path)

	def browseArchive(self):

		self.archive_path = qt.QFileDialog.getExistingDirectory(self, 'Choose archive directory')

		self.arcBox.setText(self.archive_path)

	def browseReceivers(self):

		filename = qt.QFileDialog.getOpenFileName(self, 'Open file')

		self.receiver_file = filename[0]

		self.recBox.setText(self.receiver_file)

class LocalInputDialogue(qt.QDialog):

	def __init__(self):
		super(LocalInputDialogue, self).__init__()

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/local_input_dialogue.ui', self)

		self.inpButton.clicked.connect(self.browseLocalFile)

		self.setWindowTitle('PySplit - Local catalogue input type')
		self.show()

	def browseLocalFile(self):

		filename = qt.QFileDialog.getOpenFileName(self, 'Open file')

		self.input_file = filename[0]

		self.inpBox.setText(self.input_file)

class DefaultFilterDialogue(qt.QDialog):

	def __init__(self):
		super(DefaultFilterDialogue, self).__init__()

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/default_filter_dialogue.ui', self)

		self.setWindowTitle('PySplit - Set default filter')
		self.show()

class PickingWindow(qt.QMainWindow):

	def __init__(self, catalogue):
		super(PickingWindow, self).__init__()

		# Parse args
		self.catalogue = catalogue
		self.catalogue_path = self.catalogue.path

		self.initUI()

	def initUI(self):
		uic.loadUi('ui_files/trace_window.ui', self)

		self.setWindowTitle('PySplit - Shear Wave Splitting software')
		self.show()

	def connect(self):
		# Connect each canvas to track motion - this will create a 
		self.z_comp_plot.canvas.mpl_connect('motion_notify_event', self._on_move)
		self.n_comp_plot.canvas.mpl_connect('motion_notify_event', self._on_move)
		self.e_comp_plot.canvas.mpl_connect('motion_notify_event', self._on_move)

		# Connect to clicks
		self.z_comp_plot.canvas.mpl_connect('button_press_event', self._on_click)
		self.n_comp_plot.canvas.mpl_connect('button_press_event', self._on_click)
		self.e_comp_plot.canvas.mpl_connect('button_press_event', self._on_click)

		# Connect to filter button
		self.button_applyFilter.clicked.connect(self.applyFilter)

	def _on_move(self, event):
		z_canvas = self.z_comp_plot.canvas
		n_canvas = self.n_comp_plot.canvas
		e_canvas = self.e_comp_plot.canvas
		if event.inaxes is z_canvas.ax or n_canvas.ax or e_canvas.ax:
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

			# redraw just the current cursor
			z_canvas.ax.draw_artist(self.zcursor)
			n_canvas.ax.draw_artist(self.ncursor)
			e_canvas.ax.draw_artist(self.ecursor)

			if self.w_beg_line:
				# Redraw the window line
				z_canvas.ax.draw_artist(self.z_window_beg)
				n_canvas.ax.draw_artist(self.n_window_beg)
				e_canvas.ax.draw_artist(self.e_window_beg)

			if self.pick_line:
				# Redraw the pick line
				z_canvas.ax.draw_artist(self.z_pick)
				n_canvas.ax.draw_artist(self.n_pick)
				e_canvas.ax.draw_artist(self.e_pick)

			if self.w_end_line:
				# Redraw the window line
				z_canvas.ax.draw_artist(self.z_window_end)
				n_canvas.ax.draw_artist(self.n_window_end)
				e_canvas.ax.draw_artist(self.e_window_end)				

			# blit just the redrawn area
			z_canvas.blit(z_canvas.ax.bbox)
			n_canvas.blit(n_canvas.ax.bbox)
			e_canvas.blit(e_canvas.ax.bbox)

	def _on_click(self, event):
		z_canvas = self.z_comp_plot.canvas
		n_canvas = self.n_comp_plot.canvas
		e_canvas = self.e_comp_plot.canvas		
		# Left-clicking handles the window start time
		if event.button == 1:
			# Set the window start line to be redrawn on move
			self.w_beg_line = True

			# Add the window beginning to the event stats
			self.event._add_stat("window_beg", event.xdata)

			# Make a vertical line artist
			self.z_window_beg = z_canvas.ax.axvline(event.xdata, linewidth=1, color="green", animated=True)
			self.n_window_beg = n_canvas.ax.axvline(event.xdata, linewidth=1, color="green", animated=True)
			self.e_window_beg = e_canvas.ax.axvline(event.xdata, linewidth=1, color="green", animated=True)

			# Restore the background region
			z_canvas.restore_region(self.z_background)
			n_canvas.restore_region(self.n_background)
			e_canvas.restore_region(self.e_background)

			# Redraw the window line
			z_canvas.ax.draw_artist(self.z_window_beg)
			n_canvas.ax.draw_artist(self.n_window_beg)
			e_canvas.ax.draw_artist(self.e_window_beg)

			if self.pick_line:
				# Redraw the pick line
				z_canvas.ax.draw_artist(self.z_pick)
				n_canvas.ax.draw_artist(self.n_pick)
				e_canvas.ax.draw_artist(self.e_pick)

			if self.w_end_line:
				# Redraw the window line
				z_canvas.ax.draw_artist(self.z_window_end)
				n_canvas.ax.draw_artist(self.n_window_end)
				e_canvas.ax.draw_artist(self.e_window_end)	

			# blit the redrawn area
			z_canvas.blit(z_canvas.ax.bbox)
			n_canvas.blit(n_canvas.ax.bbox)
			e_canvas.blit(e_canvas.ax.bbox)

		# Middle-clicking handles the arrival pick time
		if event.button == 2:
			# Set the pick line to be redrawn on move
			self.pick_line = True

			# Add the pick to the event stats
			self.event._add_stat("pick", event.xdata)

			# Make a vertical line artist
			self.z_pick = z_canvas.ax.axvline(event.xdata, linewidth=1, color="red", animated=True)
			self.n_pick = n_canvas.ax.axvline(event.xdata, linewidth=1, color="red", animated=True)
			self.e_pick = e_canvas.ax.axvline(event.xdata, linewidth=1, color="red", animated=True)

			# Restore the background region
			z_canvas.restore_region(self.z_background)
			n_canvas.restore_region(self.n_background)
			e_canvas.restore_region(self.e_background)

			if self.w_beg_line:
				# Redraw the window line
				z_canvas.ax.draw_artist(self.z_window_beg)
				n_canvas.ax.draw_artist(self.n_window_beg)
				e_canvas.ax.draw_artist(self.e_window_beg)

			# Redraw the pick line
			z_canvas.ax.draw_artist(self.z_pick)
			n_canvas.ax.draw_artist(self.n_pick)
			e_canvas.ax.draw_artist(self.e_pick)

			if self.w_end_line:
				# Redraw the window line
				z_canvas.ax.draw_artist(self.z_window_end)
				n_canvas.ax.draw_artist(self.n_window_end)
				e_canvas.ax.draw_artist(self.e_window_end)

			# blit the redrawn area
			z_canvas.blit(z_canvas.ax.bbox)
			n_canvas.blit(n_canvas.ax.bbox)
			e_canvas.blit(e_canvas.ax.bbox)

		# Right-clicking handles the window end time
		if event.button == 3:
			# Set the window end line to be redrawn on move
			self.w_end_line = True

			# Add the window ending to the event stats
			self.event._add_stat("window_end", event.xdata)

			# Make a vertical line artist
			self.z_window_end = z_canvas.ax.axvline(event.xdata, linewidth=1, color="green", animated=True)
			self.n_window_end = n_canvas.ax.axvline(event.xdata, linewidth=1, color="green", animated=True)
			self.e_window_end = e_canvas.ax.axvline(event.xdata, linewidth=1, color="green", animated=True)

			# Restore the background region
			z_canvas.restore_region(self.z_background)
			n_canvas.restore_region(self.n_background)
			e_canvas.restore_region(self.e_background)

			if self.w_beg_line:
				# Redraw the window line
				z_canvas.ax.draw_artist(self.z_window_beg)
				n_canvas.ax.draw_artist(self.n_window_beg)
				e_canvas.ax.draw_artist(self.e_window_beg)

			if self.pick_line:
				# Redraw the pick line
				z_canvas.ax.draw_artist(self.z_pick)
				n_canvas.ax.draw_artist(self.n_pick)
				e_canvas.ax.draw_artist(self.e_pick)

			# Redraw the window line
			z_canvas.ax.draw_artist(self.z_window_end)
			n_canvas.ax.draw_artist(self.n_window_end)
			e_canvas.ax.draw_artist(self.e_window_end)

			# blit the redrawn area
			z_canvas.blit(z_canvas.ax.bbox)
			n_canvas.blit(n_canvas.ax.bbox)
			e_canvas.blit(e_canvas.ax.bbox)

	def applyFilter(self):
		min_freq = float(self.input_minFreq.text())
		max_freq = float(self.input_maxFreq.text())

		filt = [min_freq, max_freq]

		self.plot_event(self.event_id, self.station, filt)

	def plot_event(self, event_id, station, filt=None):

		self.event_id = event_id
		self.station  = station

		# Clear the canvases
		self.z_comp_plot.canvas.ax.clear()
		self.n_comp_plot.canvas.ax.clear()
		self.e_comp_plot.canvas.ax.clear()

		# Make the window and pick lines not show on move
		self.w_beg_line = False
		self.pick_line  = False
		self.w_end_line = False

		self.z_background = None
		self.n_background = None
		self.e_background = None

		# Create an instance of the Event class
		self.event = evt.Event("{}/data/{}/event.{}.{}.*".format(self.catalogue_path, station.upper(), event_id, station.upper()))

		if not filt == None:
			self.event.filter_obspy("bandpass", filt[0], filt[1])

		self.event.plot_traces(self.z_comp_plot.canvas.ax, self.n_comp_plot.canvas.ax, self.e_comp_plot.canvas.ax)

		# Connect to each trace to grab the background once Qt has done resizing
		self.z_comp_plot.canvas.mpl_connect('draw_event', self._z_draw_event)
		self.n_comp_plot.canvas.mpl_connect('draw_event', self._n_draw_event)
		self.e_comp_plot.canvas.mpl_connect('draw_event', self._e_draw_event)

		# Draw when Qt has done all resizing
		self.z_comp_plot.canvas.draw_idle()
		self.n_comp_plot.canvas.draw_idle()
		self.e_comp_plot.canvas.draw_idle()

		# Initialise the cursor to track mouse position on the axes
		self.zcursor = self.z_comp_plot.canvas.ax.axvline(5, linewidth=1, color='0.5', animated=True)
		_, self.zydat = self.zcursor.get_data()
		self.ncursor = self.n_comp_plot.canvas.ax.axvline(5, linewidth=1, color='0.5', animated=True)
		_, self.nydat = self.ncursor.get_data()
		self.ecursor = self.e_comp_plot.canvas.ax.axvline(5, linewidth=1, color='0.5', animated=True)
		_, self.eydat = self.ecursor.get_data()

		self.connect()

		# Update the station and event information
		self._updateStationInformation(station)
		self._updateEventInformation(event_id)

	def _z_draw_event(self, evt):
		self.z_background = evt.canvas.copy_from_bbox(evt.canvas.ax.bbox)

	def _n_draw_event(self, evt):
		self.n_background = evt.canvas.copy_from_bbox(evt.canvas.ax.bbox)

	def _e_draw_event(self, evt):
		self.e_background = evt.canvas.copy_from_bbox(evt.canvas.ax.bbox)

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
		self.label_eventMag.setText(event_info.evmag.values[0])
		self.label_eventID.setText(str(event_info.sourceid.values[0]))


if __name__ == "__main__":

	app = qt.QApplication(sys.argv)
	pysplit = PySplit()
	sys.exit(app.exec_())