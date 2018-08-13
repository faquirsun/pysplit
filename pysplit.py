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

		# Specify input options and disable them on load
		self.input_options = [self.startDate_input, self.endDate_input, self.minMag_input, self.centreLon_input, self.centreLat_input, self.minRad_input, self.maxRad_input, self.tele_load_button]

		self.setWindowTitle('PySplit - Shear Wave Splitting software')
		self.setWindowIcon(QtGui.QIcon("misc/icon.png"))
		self.show()

	def connect(self):
		# File menu connections
		self.actionNew_catalogue.triggered.connect(self.newCatalogue)
		self.actionLoad_catalogue.triggered.connect(self.loadCatalogue)

		# Map connections
		self.mpl.canvas.mpl_connect('pick_event', self._on_pick)

		# Station list connection
		self.station_list.doubleClicked.connect(self.stationSelect)

		# Event list connection
		#self.event_list.doubleClicked.connect(self.pickEvent)

		# Catalogue generation button connection
		self.tele_load_button.clicked.connect(self.generate_teleseismic_catalogue)

		# Catalogue load button connection
		self.load_events_button.clicked.connect(self.load_events)
		self.load_stats_button.clicked.connect(self.load_stations)
		self.load_arrivals_button.clicked.connect(self.load_arrivals)
		self.load_waveforms_button.clicked.connect(self.load_waveforms)

		self.picked_plot_tickbox.stateChanged.connect(self.pick_state_changed)
		self.unpicked_plot_tickbox.stateChanged.connect(self.unpicked_state_changed)
		self.station_plot_tickbox.stateChanged.connect(self.station_state_changed)

	def unpicked_state_changed(self, int):
		if self.unpicked_plot_tickbox.isChecked():
			pass
			#self.catalogue.station_points.set_visible(True)
		else:
			pass
			#self.catalogue.station_points.set_visible(False)

	def station_state_changed(self, int):
		if self.station_plot_tickbox.isChecked():
			self.catalogue.station_points.set_visible(True)
		else:
			self.catalogue.station_points.set_visible(False)

	def pick_state_changed(self, int):
		if self.picked_plot_tickbox.isChecked():
			self.catalogue.pick_points.set_visible(True)
		else:
			self.catalogue.pick_points.set_visible(False)

	def _on_pick(self, event):
	    artist = event.artist
	    xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
	    xy = artist.get_offsets()
	    label = artist.get_label()
	    ind = event.ind
	    lon, lat = self.catalogue.m(xy[ind[0]][0], xy[ind[0]][1], inverse=True)

	    # If station
	    if "STAT" in label:
	    	stat = label.split(": ")[1]
	    	self._update_station_information(stat)

	    # If event
	    if "EVENT" in label:
	    	evt = label.split(": ")[1]
	    	self._update_event_information(evt)

	    	# Open an instance of Picker (TESTING)
	    	self.pickerWindow = PickingWindow(self.catalogue)
	    	self.pickerWindow.plot_event(1261, "ASK")

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
				# Set catalogue type and data source
				self.catalogue_type = "teleseismic"
				self.data_source = "IRIS"
				self.label_catTypeDisp.setText("Teleseismic")

				self.status_page.hide()
				self.teleseismic_input_page.show()
				self.station_info_page.hide()
				self.catalogue_info_page.show()

				# Create an instance of the Teleseismic Catalogue class
				self.catalogue = cat.TeleseismicCatalogue("{}/{}".format(self.catalogue_path, self.catalogue_name), self.archive_path, self.receiver_file)

			elif self.newCatalogueDialogue.local_radio.isChecked():
				# Set catalogue type
				self.catalogue_type = "local"
				self.label_catTypeDisp.setText("Local")

				self.teleseismic_input_page.hide()
				self.status_page.show()
				self.catalogue_info_page.hide()
				self.station_info_page.show()

				# Open input file dialogue
				self.localInputDialogue = LocalInputDialogue()

				if self.localInputDialogue.exec_():
					# Read in input file
					self.local_input_file = self.localInputDialogue.inpBox.text()

					# Set data source
					self.data_source = self.local_input_file

					# Create an instance of the Local Catalogue class
					self.catalogue = cat.LocalCatalogue("{}/{}".format(self.catalogue_path, self.catalogue_name), self.archive_path, self.receiver_file)

					# Create the catalogue metafile - this contains general information about the catalogue
					# and is read in when an existing catalogue is loaded
					self._generate_catalogue_metafile()
					self.generate_catalogue()

					# Plot the map of the catalogue
					self.plotter()

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


	def plotter(self):
		self.statusbar.showMessage("Plotting catalogue map...")
		self.mpl.show()
		self.mpl.canvas.ax.clear()
		self.catalogue.plot_geographic(self.mpl.canvas.ax)
		if self.station_plot_tickbox.isChecked():
			self.catalogue.plot_stations(self.mpl.canvas.ax)
		self.mpl.canvas.draw()
		self.statusbar.showMessage("Catalogue load complete.")

	def stationSelect(self, index):
		# Parse station selected
		station = self.station_list.model().data(index)

		# Load station information and print to display
		self._update_station_information(station)

		# Filter and load the corresponding events
		station_path = "{}/{}/data/{}".format(self.catalogue_path, self.catalogue_name, station.upper())

		events = glob.glob('{}/*.z'.format(station_path))

		# Populate the event list
		model = QtGui.QStandardItemModel(self.event_list)
		for event in events:
			head, tail = os.path.split(event)
			item = QtGui.QStandardItem("Event {}".format(tail.split(".")[1]))
			model.appendRow(item)
		self.event_list.setModel(model)
		self.event_list.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)

	def _update_station_information(self, station):
		station_info = self.catalogue.receiver_df.query('stat == @station')
		self.label_statName.setText(station)
		self.label_statLon.setText(f"{station_info.lon.values[0]:.4f}")
		self.label_statLat.setText(f"{station_info.lat.values[0]:.4f}")
		self.label_statElev.setText(f"{station_info.dep.values[0]:.4f}")
		self.label_statDep.setText(station_info.st_dep.values[0].isoformat().split("T")[0])
		self.label_statRet.setText(station_info.et_dep.values[0].isoformat().split("T")[0])

	def _update_event_information(self, event):
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
		self.plotter()

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
			self.plotter()

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

	def _on_move(self, event):
		if event.inaxes is self.z_comp_plot.canvas.ax:
			x = event.xdata
			self.cursor.set_data([x, x], self.ydat)
			self.cursor.set_visible(True)
			self.z_comp_plot.canvas.draw()

		if event.inaxes is self.n_comp_plot.canvas.ax:
			x = event.xdata
			self.cursor.set_data([x, x], self.ydat)
			self.cursor.set_visible(True)
			self.n_comp_plot.canvas.draw()

		if event.inaxes is self.e_comp_plot.canvas.ax:
			x = event.xdata
			self.cursor.set_data([x, x], self.ydat)
			self.cursor.set_visible(True)
			self.e_comp_plot.canvas.draw()

	def plot_event(self, event_id, station, filt=None):

		# Clear the canvases
		self.z_comp_plot.canvas.ax.clear()
		self.n_comp_plot.canvas.ax.clear()
		self.e_comp_plot.canvas.ax.clear()

		# Initialise the cursor to track mouse position on the axes
		self.cursor = self.z_comp_plot.canvas.ax.axvline(5, linewidth=1, color='0.5', visible=False)
		_, self.ydat = self.cursor.get_data()

		self.connect()

		# Update the station and event information
		#self._update_station_information(station)
		#self._update_event_information(event_id)

		# Create an instance of the Event class
		event = evt.Event("{}/data/{}/event.{}.{}.*".format(self.catalogue_path, station.upper(), event_id, station.upper()))

		if not filt == None:
			event.filter("bandpass", freqmin=filt[0], freqmax=filt[1])

		event.plot_traces(self.z_comp_plot.canvas.ax, self.n_comp_plot.canvas.ax, self.e_comp_plot.canvas.ax)

		self.z_comp_plot.canvas.draw()
		self.n_comp_plot.canvas.draw()
		self.e_comp_plot.canvas.draw()

	def _update_station_information(self, station):
		station_info = self.catalogue.receiver_df.query('stat == @station')
		self.label_statName.setText(station)
		self.label_statLon.setText(f"{station_info.lon.values[0]:.4f}")
		self.label_statLat.setText(f"{station_info.lat.values[0]:.4f}")
		self.label_statElev.setText(f"{station_info.dep.values[0]:.4f}")
		self.label_statDep.setText(station_info.st_dep.values[0].isoformat().split("T")[0])
		self.label_statRet.setText(station_info.et_dep.values[0].isoformat().split("T")[0])

	def _update_event_information(self, event):
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