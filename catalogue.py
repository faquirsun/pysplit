#!/usr/bin/env python3

"""
The Catalogue class - generate an instance of Catalogue and load existing catalogues
or create new ones from hyp files or the IRIS catalogue.

Methods
Load: loads a pre-existing catalogue, located at *path*
new_teleseismic_cat
new_local_cat
merge ???

Attributesvt
path: location of the catalogue

Author: Hemmelig
"""

# ----- Import dependencies -----
from abc import ABC, abstractmethod
import os
import sys
import glob
import ast
import pandas as pd
import datetime
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees, gps2dist_azimuth
from obspy.core import read, AttribDict
from string import Template
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

class Catalogue(ABC):
	"""
	Catalogue objects should be initiated with a path to the chosen directory.
	
	Methods
	Load catalogue (load a pre-existing catalogue)
	Teleseismic catalogue (create a new teleseismic catalogue - based on functionality from obspy)
	Local catalogue (create a new local catalogue - can be generated from hyp file)
	"""

	cmpaz  = {'N':0, 'Z':0, 'E':90}
	cmpinc = {'N':90, 'Z':0, 'E':90}

	def __init__(self):
		"""
		An abstract base class for Catalogue. The classes LocalCatalogue and TeleseismicCatalogue 
		inherit from this class.
		"""

		super().__init__()

		# Sources
		self.source_cols = ['otime', 'evlat', 'evlon', 'evdep', 'evmag', 'sourceid']
		self.source_df   = pd.DataFrame(columns=self.source_cols)

		# Arrivals
		self.arrival_cols = ['sourceid', 'receiverid', 'traveltime', 'waveform?']
		self.arrival_df   = pd.DataFrame(columns=self.arrival_cols)

		# Receivers
		self.receiver_cols = ['network', 'receiver_name', 'lat', 'lon', 'dep', 'st_dep', 'et_dep', 'receiverid']
		self.receiver_df   = pd.DataFrame(columns=self.receiver_cols)

	# --------------
	# Shared methods
	# --------------
	def newCatalogue(self, catalogue_parameters):
		self.initCatalogue(**catalogue_parameters)

		# Generate catalogue metafile
		self._generateCatalogueMetafile()

	def loadCatalogue(self, catalogue_parameters):
		self.initCatalogue(**catalogue_parameters)

	def getWaveforms(self):
		# Get window range
		window_range = self.windowRange()

		for i, arrival in self.arrival_df.iterrows():

			# If the waveforms has already been retrieved, skip it.
			if arrival['waveform?']:
				continue

			rec = self.receiver_df[self.receiver_df["receiverid"] == arrival.receiverid].receiver_name.values[0]

			for phase, ttimes in arrival["traveltime"].items():
				for traveltime in ttimes:
					# Retrieve information on the source
					otime, evlat, evlon, evdep = self._getSource(self.source_df, arrival.sourceid)

					# Calculate reference time for window to download
					window_beg = otime + float(traveltime) - window_range
					window_end = otime + float(traveltime) + window_range

					# Data existence tracker
					file_exists = True

					# Download and trim each file from the archive
					for comp in ['Z', 'N', 'E']:

						# Skip if the file doesn't exist
						if file_exists:
							pass
						else:
							continue

						evyear = str(window_beg.year)
						evjday = str(window_beg.julday)
						nxjday = str((window_beg + 86400).julday)

						# Julian day must be a 3 character string
						while (len(evjday) != 3):
							evjday = '0' + evjday
						while (len(nxjday) != 3):
							nxjday = '0' + nxjday

						# Copy all files for a given day + the next
						file_to_grab = '{}/{}/{}/*_{}_{}2.m'.format(self.archive, evyear, evjday, rec.upper(), comp)
						next_to_grab = '{}/{}/{}/*_{}_{}2.m'.format(self.archive, evyear, nxjday, rec.upper(), comp)

						# Grab files
						os.system('scp {} {}/tmp/.'.format(file_to_grab, self.catalogue_path))
						if (window_beg.julday != window_end.julday):
							os.system('scp {} {}/tmp/.'.format(next_to_grab, self.catalogue_path))

						files = sorted(glob.glob('{}/tmp/*_{}_{}2.m'.format(self.catalogue_path, rec.upper(), comp)))

						try:
							file = files[0]
						except IndexError:
							file_exists = False
							continue

						st = read('{}/tmp/*_{}_{}2.m'.format(self.catalogue_path, rec.upper(), comp))

						# Trim and save the file locally
						try:
							st.trim(window_beg, window_end)
							tr = st[0]

							# Generate name of output file
							if not os.path.exists('{}/data/{}'.format(self.catalogue_path, rec.upper())):
								os.makedirs('{}/data/{}'.format(self.catalogue_path, rec.upper()))
								os.makedirs('{}/picks/{}'.format(self.catalogue_path, rec.upper()))

							name = '{}/data/{}/source.{}.{}.{}'.format(self.catalogue_path, rec.upper(), arrival.sourceid, rec.upper(), comp.lower())

							# Write the file out to MSEED
							tr.write(name, format="MSEED")

							# Clean up variables and remove temp files
							del name
							for l in files:
								os.remove(l)
							del files

						except IndexError:
							# Clean up variables and remove temp files
							for l in files:
								os.remove(l)
							del files
							continue

					if file_exists:
						print("Data retrieved.")
						self.arrival_df.loc[i, 'waveform?'] = True

					else:
						print("The data for this arrival is not in the archive - it has been removed from the DataFrame.")
						self.arrival_df.drop(i, inplace=True)
						self.arrival_df.to_csv(self.arrival_file, index=False)

		self.arrival_df.to_csv(self.arrival_file, index=False)

	def loadSources(self):
		try:
			# Read in source file
			self.source_df = pd.read_csv(self.source_file)
			# Convert origin time strings to UTCDateTime objects
			self.source_df['otime'] = self.source_df['otime'].apply(UTCDateTime)
			return True
		except FileNotFoundError:
			return False

	def loadReceivers(self):
		try:
			# Read in receiver file
			self.receiver_df = pd.read_csv(self.receiver_file)
			# Convert start and end of deployment time strings to UTCDateTime objects
			self.receiver_df['st_dep'] = self.receiver_df['st_dep'].apply(UTCDateTime)
			self.receiver_df['et_dep'] = self.receiver_df['et_dep'].apply(UTCDateTime)
			return True		
		except FileNotFoundError:
			return False

	def loadArrivals(self):
		try:
			# Read in arrival file
			self.arrival_df = pd.read_csv(self.arrival_file)
			# Convert traveltime strings to dictionary objects
			self.arrival_df["traveltime"] = self.arrival_df["traveltime"].apply(ast.literal_eval)
			return True
		except FileNotFoundError:
		 	return False

	def loadWaveforms(self):
		# Returns true if all waveforms have been downloaded
		return self.arrival_df["waveform?"].all()

	def filterSpatial(self):
		pass

	def filterTemporal(self):
		pass

	def _getReceiver(self, df, id):
		rec  = df.query('receiverid == @id')['receiver_name'].iloc[0]
		slat = df.query('receiverid == @id')['lat'].iloc[0]
		slon = df.query('receiverid == @id')['lon'].iloc[0]
		sdep = df.query('receiverid == @id')['dep'].iloc[0]
		return rec, slat, slon, sdep

	def _getSource(self, df, id):
		otime = df.query('sourceid == @id')['otime'].iloc[0]
		evlat = df.query('sourceid == @id')['evlat'].iloc[0]
		evlon = df.query('sourceid == @id')['evlon'].iloc[0]
		evdep = df.query('sourceid == @id')['evdep'].iloc[0]
		return UTCDateTime(otime), evlat, evlon, evdep

	# --------------

	# ----------------
	# Abstract methods
	# ----------------

	@abstractmethod
	def initCatalogue(self):
		pass

	@abstractmethod
	def _generateCatalogueMetafile(self):
		pass

	@abstractmethod
	def generateCatalogue(self):
		pass

	@abstractmethod
	def getArrivals(self):
		pass

	@abstractmethod
	def windowRange(self):
		pass
	@abstractmethod
	def plotGeographic(self):
		pass

	# ----------------

class LocalCatalogue(Catalogue):
	"""

	"""

	# ----------------------
	# Class-specific methods
	# ----------------------

	def initCatalogue(self, **kwargs):
		# Parse keyword arguments
		self.catalogue_name = kwargs.get("catalogue_name")
		self.catalogue_type = kwargs.get("catalogue_type")
		self.catalogue_path = kwargs.get("catalogue_path")
		self.data_source    = kwargs.get("data_source")
		self.creation_date  = kwargs.get("creation_date")
		self.archive_path   = kwargs.get("archive_path")
		self.archive_format = kwargs.get("archive_format")
		self.receiver_file  = kwargs.get("receiver_file")
		self.start_date     = kwargs.get("start_date")
		self.end_date       = kwargs.get("end_date")
		self.local_input    = kwargs.get("local_input")

		# Check if the path already exists - otherwise, create it.
		if not os.path.exists(self.catalogue_path):
			os.makedirs(self.catalogue_path)
			os.makedirs("{}/metafiles".format(self.catalogue_path))
			os.makedirs("{}/data".format(self.catalogue_path))
			os.makedirs("{}/tmp".format(self.catalogue_path))
			os.makedirs("{}/plots".format(self.catalogue_path))
			os.makedirs("{}/picks".format(self.catalogue_path))
			print("Catalogue located at {}.".format(self.catalogue_path))

		receiver_file_new = "{}/metafiles/receivers.txt".format(self.catalogue_path)
		# Copy receiver file into metafiles directory
		os.system("cp {} {}".format(self.receiver_file, receiver_file_new))

		#self.receiver_file = 

		# Set source file location
		self.source_file = "{}/metafiles/sources.txt".format(self.catalogue_path)

		# Set arrival file location
		self.arrival_file = "{}/metafiles/arrivals.txt".format(self.catalogue_path)

		# Set receiver file location
		self.receiver_file = "{}/metafiles/receivers.txt".format(self.catalogue_path)
		self.receiver_df = pd.read_csv(self.receiver_file, sep=',')

	def _generateCatalogueMetafile(self):
		# Read in the template file
		filein = open('templates/local_metafile_template.txt')
		src = Template(filein.read())

		# Set the general metadata parameters
		d_cat = {'catalogue_name': self.catalogue_name,
			 	 'catalogue_type': self.catalogue_type,
				 'catalogue_path': self.catalogue_path,
				 'data_source': self.data_source,
				 'creation_date': self.creation_date,
				 'archive_path': self.archive_path,
				 'archive_format': self.archive_format,
				 'receiver_file': self.receiver_file,
				 'start_date': self.start_date,
				 'end_date': self.end_date,
				 'local_input': self.local_input}

		# Write out the dictionaries of the parameters - safe substitute ignores missing template variables
		# and just writes a blank
		output = src.safe_substitute(d_cat)

		# Write out the catalogue metafile
		with open('{}/metafiles/catalogue_metafile.txt'.format(self.catalogue_path), 'w') as o:
			o.write(output)

	def generateCatalogue(self):
		"""
		Extracts earthquake source information from an input file and
		creates a catalogue.

		Supported input types:
			SeisLoc directory
			PySplit catalogue
			.hyp file
			IRIS catalogue## not yet supported
		"""

		# Check that the input file/path provided exists
		if not os.path.exists(self.local_input):
			print("File doesn't exist.")
			print("Please provide a valid path.")

		# Input type: SeisLoc
		if self.data_source == "SeisLoc":
			# Parse SeisLoc directory
			sources = glob.glob("{}/outputs/*.event".format(self.local_input))

			# Create DataFrame for source information. Space preallocated
			# to improve efficiency.
			self.source_df = pd.DataFrame(index=np.arange(0, len(source_f)), columns=self.source_cols)

			# Loop through the sources and create the source DataFrame, as well
			# as the arrivals DataFrame
			for idx, source in enumerate(sources):
				# Handle source DataFrame 
				lines = []
				with open(source, 'r') as f:
					for line in f:
						lines.append(line)

				s_info = lines[1]
				s_info = s_info.rstrip().split(",")
				otime = UTCDateTime(s_info[0].replace(" ", "T"))
				evlon = s_info[2]
				evlat = s_info[3]
				evdep = s_info[4]
				evmag = "?"

				self.source_df.loc[idx] = [otime, evlat, evlon, evdep, evmag, idx]

				# Handle arrival DataFrame
				arrival = glob.glob("{}/outputs/{}.stn".format(self.local_input, source[:-6]))

				lines = []
				with open(arrival, 'r') as f:
					for line in f:
						lines.append(line)
				lines = lines[1:]

				for i, line in enumerate(lines):
					line.rstrip.split(",")

					receiver = line[0]
					phase   = line[1]
					model_t = UTCDateTime(line[2].replace(" ", "T"))
					if line[3] == -1.0:
						pass
					else:
						pick_t = UTCDateTime(line[3].replace(" ", "T"))
					if line[4] == -1.0:
						pick_e = None

					ttime = model_t - otime

					self.arrival_df.loc[i + 1 * idx] = [idx, self._lookupReceiverID(receiver), ttime, False]

		# Input type: PySplit catalogues
		elif self.data_source == "PySplit":
			# Copy the input file across into the catalogue directory
			head, tail = os.path.split(self.local_input)
			input_file = "{}/{}/metafiles/{}".format(self.catalogue_path, self.catalogue_name, tail)
			os.system("cp {} {}".format(self.local_input, input_file))

			self.source_df = pd.read_csv(input_file)
			self.source_df["otime"] = self.source_df["otime"].apply(UTCDateTime)
			self.source_df["sourceid"] = self.source_df.index

		# Input type: .hyp file
		elif self.data_source == ".hyp file":
			# Copy the input file across into the catalogue directory
			head, tail = os.path.split(self.local_input)
			input_file = "{}/metafiles/{}".format(self.catalogue_path, tail)
			os.system("cp {} {}".format(self.local_input, input_file))

			if tail[-3:] == "hyp":
				lines = []
				with open(input_file, 'r') as f:
					for line in f:
						if "GEOGRAPHIC" in line:
							lines.append(line)

				# Create DataFrame for source information. Space preallocated to
				# improve efficiency.
				self.source_df = pd.DataFrame(index=np.arange(0, len(lines)), columns=self.source_cols)

				for i in range(len(lines)):
					line = lines[i].rstrip().split()

					# Get origin time
					evyear, evmonth, evday  = int(line[2]), int(line[3]), int(line[4])
					evhour, evmin, evsecflt = int(line[5]), int(line[6]), float(line[7])
					evsec  = int(evsecflt)
					evmsec = int(1e6 * (evsecflt - evsec))
					otime  = UTCDateTime(evyear, evmonth, evday, evhour, evmin, evsec, evmsec)

					# Get hypocentral location
					evlat = float(line[9])
					evlon = float(line[11])
					evdep = float(line[13])
					evmag = "?" #### ADD PROPER HANDLING FOR HYP FILES THAT HAVE MAGNITUDES IN

					self.source_df.loc[i] = [otime, evlat, evlon, evdep, evmag, i]

			else:
				print("Please provide a .hyp file.")

		# Input type: IRIS catalogue format
		elif self.data_source == "IRIS":
			# Add support here
			# 
			#
			pass

		# Output the catalogue
		self.source_df.to_csv(self.source_file, index=False)

	def getArrivals(self, phases=None, input_file=None, input_type=None):
		"""
		Extracts seismic phase arrival times from an input file and 
		creates a list of arrivals.
		"""
		if input_file == None:
			print("You've not provided file.")
			return

		# Input type: PySplit
		if input_type == "PySplit":
			# Read in the input file
			sources_df = pd.read_csv(input_file)
			# The arrival file simply needs to contain the origin time
			# The waveform script will give a +/- 60 second window
			# around this
			self.arrival_df = pd.DataFrame(index=np.arange(0, len(self.receiver_df) * len(sources_df)), columns=self.arrival_cols)
			for idx, source in sources_df.iterrows():
				for jdx, receiver in self.receiver_df.iterrows():
					# Want to check if the receiver was available at the time of the source
					otime = UTCDateTime(source["otime"])
					stdp  = UTCDateTime(receiver["st_dep"])
					etdp  = UTCDateTime(receiver["et_dep"])
					if (otime >= stdp) & (otime <= etdp):
						self.arrival_df.loc[jdx + idx * 1] = [idx, receiver['receiverid'], "{'arr': [0.0]}", False]
					else:
						self.arrival_df.loc[jdx + idx * 1] = ["-", "-", "-", "-"]
						continue

		# Input type: .hyp
		if input_type == ".hyp":
			lines = []
			S_lines = []
			with open(input_file, 'r') as f:
				for line in f:
					if ">" in line:
						if " P " in line:
							continue
						else:
							lines.append(line)

					if (">" in line) and (" S " in line):
						S_lines.append(line)

			# Preallocate space for the arrivals DataFrame
			self.arrival_df = pd.DataFrame(index=np.arange(0, len(S_lines)), columns=self.arrival_cols)

			sourceid = -1
			# Every time a line starting with PHASE is passed, increment the sourceid
			for i in range(len(lines)):
				line = lines[i].rstrip().split()

				if line[0] == "PHASE":
					sourceid += 1
					continue

				else:
					idx = i - (sourceid + 1)
					print(line[0])
					receiver = line[0]
					try:
						self.arrival_df.loc[idx] = [sourceid, self._lookupReceiverID(receiver), "{'arr': [line[15]]}", False]
					except:
						print(receiver, " not found, check it exists?")
						self.arrival_df.loc[idx] = ["-", "-", "-", "-"]
						continue

		# Remove all sources that didn't have complete information
		self.arrival_df = self.arrival_df.drop(self.arrival_df[self.arrival_df.sourceid == "-"].index)

		# Reset the index
		self.arrival_df.index = pd.RangeIndex(len(self.arrival_df.index))

		# Output the catalogue
		self.arrival_df.to_csv(self.arrival_file, index=False)

	def windowRange(self):
		return 60.

	def plotGeographic(self, map_widget, lon0=None, lon1=None, lat0=None, lat1=None):

		# Find appropriate geographical region
		lons = self.source_df.evlon.values
		lats = self.source_df.evlat.values

		if lon0 == None: 
			londiff = (lons.max() - lons.min()) * 0.1
			latdiff = (lats.max() - lats.min()) * 0.1
			self.lon0 = lons.min() - londiff
			self.lon1 = lons.max() + londiff
			self.lat0 = lats.min() - latdiff
			self.lat1 = lats.max() + latdiff

		else:
			self.lon0 = lon0
			self.lon1 = lon1
			self.lat0 = lat0
			self.lat1 = lat1
			
		self.map = map_widget.canvas._localMap(lon0=self.lon0, lat0=self.lat0, lon1=self.lon1, lat1=self.lat1)

		# Set plot details (axes labels etc)
		self.map.set_xlabel("Longitude, degrees", fontsize=10)
		self.map.set_ylabel("Latitude, degrees", fontsize=10)
		self.map.set_extent([self.lon0, self.lon1, self.lat0, self.lat1], ccrs.PlateCarree())

		# Try rescaling the image now
		self.map.set_aspect('auto')

		tolerance = 10
		for i in range(len(lons)):
			self.map.scatter(lons[i], lats[i], 12, marker='o', color='k', picker=tolerance, zorder=10, label="SOURCE: {}".format(self.source_df.sourceid[i]))

	def plotReceivers(self, map_widget):
		lons = self.receiver_df.lon.values
		lats = self.receiver_df.lat.values
		tolerance = 10
		for i in range(len(lons)):
			self.map.scatter(lons[i], lats[i], 75, marker="v", color="green", picker=tolerance, zorder=15, label="REC: {}".format(self.receiver_df.receiver_name[i]))

	def _lookupReceiverID(self, receiver):
		# Returns the receiver ID for given receiver name
		return self.receiver_df.query('receiver_name == @receiver')["receiverid"].iloc[0]



class TeleseismicCatalogue(Catalogue):
	"""

	"""

	def initCatalogue(self, **kwargs):
		# Parse keyword arguments
		self.catalogue_name = kwargs.get("catalogue_name")
		self.catalogue_type = kwargs.get("catalogue_type")
		self.catalogue_path = kwargs.get("catalogue_path")
		self.data_source    = kwargs.get("data_source")
		self.creation_date  = kwargs.get("creation_date")
		self.archive_path   = kwargs.get("archive_path")
		self.archive_format = kwargs.get("archive_format")
		self.receiver_file  = kwargs.get("receiver_file")
		self.start_date     = kwargs.get("start_date")
		self.end_date       = kwargs.get("end_date")
		self.local_input    = kwargs.get("local_input")
		self.minmag         = kwargs.get("minmag")
		self.clon           = kwargs.get("clon")
		self.clat           = kwargs.get("clat")
		self.minrad         = kwargs.get("minrad")
		self.maxrad         = kwargs.get("maxrad")


		# Check if the path already exists - otherwise, create it.
		if not os.path.exists(self.catalogue_path):
			os.makedirs(self.catalogue_path)
			os.makedirs("{}/metafiles".format(self.catalogue_path))
			os.makedirs("{}/data".format(self.catalogue_path))
			os.makedirs("{}/tmp".format(self.catalogue_path))
			os.makedirs("{}/plots".format(self.catalogue_path))
			os.makedirs("{}/picks".format(self.catalogue_path))
			print("Catalogue located at {}.".format(self.catalogue_path))

		# Copy receiver file into metafiles directory
		os.system("cp {} {}/metafiles/receivers.txt".format(self.receiver_file, self.catalogue_path))

		# Set source file location
		self.source_file = "{}/metafiles/sources.txt".format(self.catalogue_path)

		# Set arrival file location
		self.arrival_file = "{}/metafiles/arrivals.txt".format(self.catalogue_path)

		# Read in receiver file
		self.receiver_file = "{}/metafiles/receivers.txt".format(self.catalogue_path)
		self.receiver_df = pd.read_csv(self.receiver_file, sep=',')

	def _generateCatalogueMetafile(self):
		# Read in the template file
		filein = open('templates/teleseismic_metafile_template.txt')
		src = Template(filein.read())

		# Set the general metadata parameters
		d_cat = {'catalogue_name': self.catalogue_name,
			 	 'catalogue_type': self.catalogue_type,
				 'catalogue_path': self.catalogue_path,
				 'data_source': self.data_source,
				 'creation_date': self.creation_date,
				 'archive_path': self.archive_path,
				 'archive_format': self.archive_format,
				 'receiver_file': self.receiver_file,
				 'start_date': self.start_date,
				 'end_date': self.end_date,
				 'minmag': self.minmag,
				 'clon': self.clon,
				 'clat': self.clat,
				 'minrad': self.minrad,
				 'maxrad': self.maxrad}

		# Write out the dictionaries of the parameters - safe substitute ignores missing template variables
		# and just writes a blank
		output = src.safe_substitute(d_cat)

		# Write out the catalogue metafile
		with open('{}/metafiles/catalogue_metafile.txt'.format(self.catalogue_path), 'w') as o:
			o.write(output)

	def generateCatalogue(self):
		# Initialise a client object
		client = Client("IRIS")

		# Convert provided start and end times to UTCDateTime
		st = UTCDateTime(self.start_date)
		et = UTCDateTime(self.end_date)

		sources = client.get_events(starttime=st, endtime=et, minmagnitude=self.minmag,
									latitude=self.clat, longitude=self.clon, 
									minradius=self.minrad, maxradius=self.maxrad)

		# Create DataFrame for source information. Space preallocated to
		# improve efficiency.
		self.source_df = pd.DataFrame(index=np.arange(0, sources.count()), columns=self.source_cols)

		for i in range(sources.count()):
			try:
				source = sources[i]
				otime = source.preferred_origin().get('time')
				evlat = source.preferred_origin().get('latitude')
				evlon = source.preferred_origin().get('longitude')
				evdep = source.preferred_origin().get('depth') / 1000.0
				evmag = source.preferred_magnitude().get('mag')

				self.source_df.loc[i] = [otime, evlat, evlon, evdep, evmag, i]

			except TypeError:
				print("No recorded depth for source:")
				print("     Origin time: {}".format(otime))
				print("        Latitude: {}".format(evlat))
				print("       Longitude: {}".format(evlon))
				print("Removed from catalogue.")

				self.source_df.loc[i] = ["-", "-", "-", "-", "-", i]

		# Remove all sources that didn't have complete information
		self.source_df = self.source_df.drop(self.source_df[self.source_df.otime == "-"].index)

		# Reset the index
		self.source_df.index = pd.RangeIndex(len(self.source_df.index))

		# Output the catalogue
		self.source_df.to_csv(self.source_file, index=False)

	def getArrivals(self, phases=["SKS"], input_file=None):
		# Initialise the velocity model
		model = TauPyModel(model="ak135")

		# Evaluate the arrivals
		for i, receiver in self.receiver_df.iterrows():
			print("i", i)
			print(" ")
			print(receiver)
			# Filter for sources occurring during receiver deployment
			tmp_df = self.source_df[self.source_df['otime'].between(receiver.st_dep, receiver.et_dep)]

			for j, source in tmp_df.iterrows():
				print("j", j)
				tmp_dist = locations2degrees(source.evlat, source.evlon, receiver.lat, receiver.lon)
				phase_arr = model.get_travel_times(source_depth_in_km=source.evdep,
												   distance_in_degree=tmp_dist,
												   phase_list=phases,
												   receiver_depth_in_km=receiver.dep)

				# Check if any of the requested phase arrivals exist
				if not phase_arr:
					print("No arrivals found for the following phases:")
					for i in range(len(phases)):
						print("     {}".format(phases[i]))
					continue

				else:
					# Create dictionary for phase traveltimes
					phase_dict = {}
					for i in range(len(phase_arr)):
						key  = phase_arr[i].purist_name
						time = phase_arr[i].time
						
						print("Adding {} arrival to dictionary:".format(key))
						if key in phase_dict:
							phase_dict[key].append(phase_arr[i].time)
						else:
							phase_dict[key] = []
							phase_dict[key].append(time) 

					tmp = pd.DataFrame([int(source.sourceid), int(receiver.receiverid), phase_dict, False], index=self.arrival_cols).T
					self.arrival_df = self.arrival_df.append(tmp, ignore_index=True)

		self.arrival_df.to_csv(self.arrival_file, index=False)

		# Now filter the receivers against this list
		available_receivers = self.arrival_df.drop_duplicates(subset="receiverid", keep="first").receiverid
		available_receivers_names = []

		for i, receiverid in available_receivers.items():
			rec, slat, slon, sdep = self._getReceiver(self.receiver_df, receiverid)
			available_receivers_names.append(rec)

		self.receiver_df = self.receiver_df.loc[self.receiver_df['receiver_name'].isin(available_receivers_names)]

		# Output the new receiver file
		self.receiver_df.to_csv(self.receiver_file, index=False)

	def windowRange(self):
		return 300.

	def plotGeographic(self, map_widget):
		# Find the centre of the array
		lon_cen = self.receiver_df.lon.mean()
		lat_cen = self.receiver_df.lat.mean()

		# Retrieve the source longitudes and latitudes
		lons = self.source_df.evlon.values
		lats = self.source_df.evlat.values

		ax = map_widget.canvas._teleseismicMap(lon_cen, lat_cen)

		#ax.add_wms(wms='http://vmap0.tiles.osgeo.org/wms/vmap0', layers=['basic'])
		tolerance = 10
		for i in range(len(lons)):
			ax.scatter(lons[i], lats[i], 12, marker='o', color='k', picker=tolerance, zorder=10, label="SOURCE: {}".format(self.source_df.sourceid[i]),  transform=ccrs.Geodetic())