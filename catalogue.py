#!/usr/bin/env python3

"""
The Catalogue class - generate an instance of Catalogue and load existing catalogues
or create new ones from hyp files or the IRIS catalogue.

Methods
Load: loads a pre-existing catalogue, located at *path*
new_teleseismic_cat
new_local_cat
merge ???

Attributes
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
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
from obspy.core import read
from string import Template
import numpy as np
import cartopy.crs as ccrs
import pathlib
import shutil
import metainfo as psm

class Catalogue(ABC):
	"""
	Catalogue class

	Abstract Base Class for LocalCatalogue and TeleseismicCatalogue. A catalogue
	contains information pertaining to 

	Attributes
	----------
	
	Methods
	-------
	Load catalogue (load a pre-existing catalogue)
	Teleseismic catalogue (create a new teleseismic catalogue - based on functionality from obspy)
	Local catalogue (create a new local catalogue - can be generated from hyp file)
	"""

	cmpaz  = {'N':0, 'Z':0, 'E':90}
	cmpinc = {'N':90, 'Z':0, 'E':90}
	src_cols = ['otime', 'lat', 'lon', 'dep', 'mag', 'sourceid']
	arr_cols = ['sourceid', 'receiver', 'traveltime', 'waveform?']

	def __init__(self, new=False, **kwargs):
		"""
		Class initialisation

		Parameters
		----------
		new : bool, optional
			If True, will also generate a metafile (default: False)
		kwargs : Dict
			Dictionary containing Catalogue attributes and values
		"""

		super().__init__()

		self.src_df = pd.DataFrame(columns=self.src_cols)
		self.arr_df = pd.DataFrame(columns=self.arr_cols)

		for key, value in kwargs.items():
			setattr(self, key, value)

		self.network = psm.Network(self.rec_file)

		self.cat_dir = pathlib.Path(self.catalogue_path)
		dirs = ["", "metafiles", "data", "tmp", "plots", "picks"]
		for dir_ in dirs:
			new_dir = self.cat_dir / dir_
			new_dir.mkdir(exist_ok=True)

		new_rec_file = self.cat_dir / "metafiles" / "receivers.txt"
		rec_file     = pathlib.Path(self.rec_file)
		try:
			shutil.copy(str(rec_file), str(new_rec_file))
		except shutil.SameFileError:
			pass

		self.src_file = self.cat_dir / "metafiles" / "sources.txt"
		self._loadSources()
		self.arr_file = self.cat_dir / "metafiles" / "arrivals.txt"
		#self.loadArrivals()

		if new:
			self._generateCatalogueMetafile()
			self._generateCatalogue()

	def getWaveforms(self):
		"""
		Retrieves waveform data from an archive

		"""
		for i, arrival in self.arr_df.iterrows():
			if arrival['waveform?']:
				continue

			rec = self.network.lookupRec(arrival.receiverid)
			src = self.lookupSrc(arrival.sourceid)
			data_dir = self.cat_dir / "data" / rec.name.upper()
			data_dir.mkdir(exist_ok=True, parents=True)
			pick_dir = self.cat_dir / "picks" / rec.name.upper()
			pick_dir.mkdir(exist_ok=True, parents=True)
			tmp_dir  = self.cat_dir / "tmp"

			for phase, ttimes in arrival["traveltime"].items():
				for traveltime in ttimes:
					wbeg = src.otime + float(traveltime) - self.window
					wend = src.otime + float(traveltime) + self.window

					evyear = str(wbeg.year)
					evjday = str(wbeg.julday).zfill(3)
					nxjday = str((wbeg + 86400).julday).zfill(3)

					evfile = self.cat_dir / self.archive_format.format(year=evyear, jday=evjday, receiver=rec.name.upper(), comp="*")
					nxfile = self.cat_dir / self.archive_format.format(year=evyear, jday=nxjday, receiver=rec.name.upper(), comp="*")

					os.system('scp {} {}/tmp/.'.format(evfile, self.catalogue_path))
					if (wbeg.julday != wend.julday):
						os.system('scp {} {}/tmp/.'.format(nxfile, self.catalogue_path))

					files = list(tmp_dir.glob("*"))

					# Only want to use data where we have all 3 components
					# The first test ensures we also capture any next day files
					if (len(files) % 3 != 0) and (len(files) != 0):
						self.arr_df.drop(i, inplace=True)
						continue

					st = obspy.Stream()
					for file in files:
						st += read(str(file))

					try:
						st.trim(wbeg, wend)
					except IndexError:
						print("Window require not within time frame of data.")
						for file in files:
							file.unlink()
						continue

					for tr in st:
						name = '{}/data/{}/source.{}.{}.{}'.format(self.catalogue_path, rec.name.upper(), 
																   arrival.sourceid, rec.upper(), comp.lower())
						tr.write(name, format="MSEED")

					for file in files:
						file.unlink()

					self.arr_df.loc[i, 'waveform?'] = True

			if (i % 10) == 0:
				self.arr_df.to_csv(self.arr_file, index=False)

		self.arr_df.to_csv(self.arr_file, index=False)

	def loadSources(self):
		if self.src_file.is_file():
			self.src_df = pd.read_csv(self.src_file)
			self.src_df['otime'] = self.src_df['otime'].apply(UTCDateTime)

	def loadArrivals(self):
		if self.arr_file.is_file():
			self.arr_df = pd.read_csv(self.arr_file)
			self.arr_df["traveltime"] = self.arr_df["traveltime"].apply(ast.literal_eval)

			return True

		else:
			return False

	def loadWaveforms(self):
		# Returns true if all waveforms have been downloaded
		return self.arr_df["waveform?"].all()

	def lookupSrc(self, sourceid):
		"""
		Queries the source DataFrame for 'sourceid'.

		Parameters
		----------
		source : str
			Source name/unique identifier

		Returns
		-------
		Source object
		"""
		if not type(sourceid) == str:
			print("You must identify the source with a string.")
			return
		try:
			src = self.src_df.loc[self.src_df["sourceid"] == sourceid]
			return Source(src)
		except IndexError:
			print("There is no source with id \"{}\" in this catalogue".format(sourceid))

	# --------------

	@property
	def window(self):
		"""Get and set the time window padding either side of the arrival of
		the signal at the receiver. Used when slicing mSEED files. (units: s)"""
		return self._window

	@window.setter
	def window(self, value):
		self._window = value
	
	# ----------------
	# Abstract methods
	# ----------------

	@abstractmethod
	def _generateCatalogueMetafile(self):
		pass

	@abstractmethod
	def _generateCatalogue(self):
		pass

	@abstractmethod
	def getArrivals(self):
		pass

	@abstractmethod
	def plotGeographic(self):
		pass

	# ----------------


class LocalCatalogue(Catalogue):
	"""LocalCatalogue class

	Inherits from Catalogue class, providing functionality specific to catalogues
	of local sources

	Attributes
	----------

	Methods
	-------

	Properties
	----------

	"""

	# ----------------------
	# Class-specific methods
	# ----------------------

	window = 60.

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
				 'archive_type': self.archive_type,
				 'archive_format': self.archive_format,
				 'rec_file': self.rec_file,
				 'start_date': self.start_date,
				 'end_date': self.end_date,
				 'local_input': self.local_input}

		# Write out the dictionaries of the parameters - safe substitute ignores missing template variables
		# and just writes a blank
		output = src.safe_substitute(d_cat)

		meta_path = pathlib.Path(self.catalogue_path) / "metafiles" / "catalogue_metafile.txt"
		with open(meta_path, 'w') as f:
			f.write(output)

	def _generateCatalogue(self):
		"""
		Extracts earthquake source information from an input file and
		creates a catalogue.

		Supported input types:
			SeisLoc directory
			PySplit catalogue
			.hyp file
			IRIS catalogue - not yet supported
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
			self.src_df = pd.DataFrame(index=np.arange(0, len(source_f)), columns=self.src_cols)

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
				lon = s_info[2]
				lat = s_info[3]
				dep = s_info[4]
				mag = "?"

				self.src_df.loc[idx] = [otime, lat, lon, dep, mag, idx]

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

					self.arr_df.loc[i + 1 * idx] = [idx, self._lookupReceiverID(receiver), ttime, False]

		# Input type: PySplit catalogues
		elif self.data_source == "PySplit":
			# Copy the input file across into the catalogue directory
			head, tail = os.path.split(self.local_input)
			input_file = "{}/{}/metafiles/{}".format(self.catalogue_path, self.catalogue_name, tail)
			os.system("cp {} {}".format(self.local_input, input_file))

			self.src_df = pd.read_csv(input_file)
			self.src_df["otime"] = self.src_df["otime"].apply(UTCDateTime)
			self.src_df["sourceid"] = self.src_df.index

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
				self.src_df = pd.DataFrame(index=np.arange(0, len(lines)), columns=self.src_cols)

				for i in range(len(lines)):
					line = lines[i].rstrip().split()

					# Get origin time
					evyear, evmonth, evday  = int(line[2]), int(line[3]), int(line[4])
					evhour, evmin, evsecflt = int(line[5]), int(line[6]), float(line[7])
					evsec  = int(evsecflt)
					evmsec = int(1e6 * (evsecflt - evsec))
					otime  = UTCDateTime(evyear, evmonth, evday, evhour, evmin, evsec, evmsec)

					# Get hypocentral location
					lat = float(line[9])
					lon = float(line[11])
					dep = float(line[13])
					mag = "?" #### ADD PROPER HANDLING FOR HYP FILES THAT HAVE MAGNITUDES IN

					self.src_df.loc[i] = [otime, lat, lon, dep, mag, i]

			else:
				print("Please provide a .hyp file.")

		# Input type: IRIS catalogue format
		elif self.data_source == "IRIS":
			# Add support here
			# 
			#
			pass

		# Output the catalogue
		self.src_df.to_csv(self.src_file, index=False)

	def getArrivals(self, phases=None, input_file=None, input_type=None):
		"""
		Extracts seismic phase arrival times from an input file and 
		creates a list of arrivals.

		Parameters
		----------
		phases : list, optional
			List of phase identification strings
		input_file : 

		input_type : , optional
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
			self.arr_df = pd.DataFrame(index=np.arange(0, len(self.network.receivers) * len(sources_df)), columns=self.arr_cols)
			for idx, source in sources_df.iterrows():
				for jdx, receiver in self.network.receivers.iterrows():
					# Want to check if the receiver was available at the time of the source
					otime = UTCDateTime(source["otime"])
					stdp  = UTCDateTime(receiver["st_dep"])
					etdp  = UTCDateTime(receiver["et_dep"])
					if (otime >= stdp) & (otime <= etdp):
						self.arr_df.loc[jdx + idx * 1] = [idx, receiver['receiverid'], "{'arr': [0.0]}", False]
					else:
						self.arr_df.loc[jdx + idx * 1] = ["-", "-", "-", "-"]
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
			self.arr_df = pd.DataFrame(index=np.arange(0, len(S_lines)), columns=self.arr_cols)

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
						self.arr_df.loc[idx] = [sourceid, self._lookupReceiverID(receiver), "{'arr': [line[15]]}", False]
					except:
						print(receiver, " not found, check it exists?")
						self.arr_df.loc[idx] = ["-", "-", "-", "-"]
						continue

		# Remove all sources that didn't have complete information
		self.arr_df = self.arr_df.drop(self.arr_df[self.arr_df.sourceid == "-"].index)
		self.arr_df.index = pd.RangeIndex(len(self.arr_df.index))
		self.arr_df.to_csv(self.arr_file, index=False)

		available_receivers = self.network.filterReceivers(self.arr_df)
		available_receivers.to_csv(self.rec_file, index=False)

	def plotGeographic(self, map_widget, lims=None):

		# Find appropriate geographical region
		lons = self.src_df.lon.values
		lats = self.src_df.lat.values
		sids = self.src_df.sourceid.values

		if lims == None: 
			londiff = (lons.max() - lons.min()) * 0.1
			latdiff = (lats.max() - lats.min()) * 0.1
			self.lon0 = lons.min() - londiff
			self.lon1 = lons.max() + londiff
			self.lat0 = lats.min() - latdiff
			self.lat1 = lats.max() + latdiff

		else:
			self.lon0 = lims["lon0"]
			self.lon1 = lims["lon1"]
			self.lat0 = lims["lat0"]
			self.lat1 = lims["lat1"]
			
		self.map = map_widget.canvas._localMap(lon0=self.lon0, lat0=self.lat0, lon1=self.lon1, lat1=self.lat1)

		# Set plot details (axes labels etc)
		self.map.set_xlabel("Longitude, degrees", fontsize=10)
		self.map.set_ylabel("Latitude, degrees", fontsize=10)
		self.map.set_extent([self.lon0, self.lon1, self.lat0, self.lat1], ccrs.PlateCarree())

		# Try rescaling the image now
		self.map.set_aspect('auto')

		tolerance = 10
		for i in range(len(lons)):
			self.map.scatter(lons[i], lats[i], 12, marker='o', color='k', picker=tolerance, zorder=10, label="SOURCE: {}".format(sids[i]))

	def plotReceivers(self, map_widget):
		lons = self.network.receivers.lon.values
		lats = self.network.receivers.lat.values
		tolerance = 10
		for i in range(len(lons)):
			self.map.scatter(lons[i], lats[i], 75, marker="v", color="green", picker=tolerance, zorder=15, label="REC: {}".format(self.network.receivers.name[i]))

	def _lookupReceiverID(self, receiver):
		# Returns the receiver ID for given receiver name
		return self.network.receivers.query('name == @receiver')["receiverid"].iloc[0]


class TeleseismicCatalogue(Catalogue):
	"""TeleseismicCatalogue class

	Inherits from Catalogue class, providing functionality specific to catalogues
	of teleseismic sources

	Attributes
	----------
	window : float
		Default 

	Methods
	-------

	Properties
	----------

	"""

	window = 300.

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
				 'archive_type': self.archive_type,
				 'archive_format': self.archive_format,
				 'rec_file': self.rec_file,
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

	def _generateCatalogue(self):
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
		self.src_df = pd.DataFrame(index=np.arange(0, sources.count()), columns=self.src_cols)

		for i in range(sources.count()):
			try:
				source = sources[i]
				otime = source.preferred_origin().get('time')
				lat = source.preferred_origin().get('latitude')
				lon = source.preferred_origin().get('longitude')
				dep = source.preferred_origin().get('depth') / 1000.0
				mag = source.preferred_magnitude().get('mag')

				# sourceid = ### Come up with file naming format
				sourceid = i

				self.src_df.loc[i] = [otime, lat, lon, dep, mag, sourceid]

			except TypeError:
				print("No recorded depth for source:")
				print("     Origin time: {}".format(otime))
				print("        Latitude: {}".format(lat))
				print("       Longitude: {}".format(lon))
				print("Removed from catalogue.")

				self.src_df.loc[i] = ["-", "-", "-", "-", "-", i]

		# Remove all sources that didn't have complete information
		self.src_df = self.src_df.drop(self.src_df[self.src_df.otime == "-"].index)

		# Reset the index
		self.src_df.index = pd.RangeIndex(len(self.src_df.index))

		# Output the catalogue
		self.src_df.to_csv(self.src_file, index=False)

	def getArrivals(self, phases=["SKS"], input_file=None):
		"""

		"""
		model = TauPyModel(model="ak135")

		for i, receiver in self.network.receivers.iterrows():
			tmp_df = self.src_df[self.src_df['otime'].between(receiver.deployment, receiver.retrieval)]

			for j, source in tmp_df.iterrows():
				tmp_dist = locations2degrees(source.lat, source.lon, receiver.lat, receiver.lon)
				phase_arr = model.get_travel_times(source_depth_in_km=source.dep,
												   distance_in_degree=tmp_dist,
												   phase_list=phases,
												   receiver_depth_in_km=receiver.elv)

				if not phase_arr:
					continue
				else:
					phase_dict = {}
					for i in range(len(phase_arr)):
						key  = phase_arr[i].purist_name
						time = phase_arr[i].time
						
						if key in phase_dict:
							phase_dict[key].append(phase_arr[i].time)
						else:
							phase_dict[key] = []
							phase_dict[key].append(time)

					tmp = pd.DataFrame([int(source.sourceid), receiver.name, phase_dict, False], index=self.arr_cols).T
					self.arr_df = self.arr_df.append(tmp, ignore_index=True)

		self.arr_df.to_csv(self.arr_file, index=False)

		available_receivers = self.network.filterReceivers(self.arr_df)
		available_receivers.to_csv(self.rec_file, index=False)

	def plotGeographic(self, map_widget=None):
		network_centre = self.network.centre
		lons = self.src_df.lon.values
		lats = self.src_df.lat.values

		if not map_widget:
			pass
		else:
			ax = map_widget.canvas._teleseismicMap(network_centre[0], network_centre[1])
		tolerance = 10
		for i in range(len(lons)):
			ax.scatter(lons[i], lats[i], 12, marker='o', color='k', picker=tolerance, zorder=10, label="SOURCE: {}".format(self.src_df.sourceid[i]),  transform=ccrs.Geodetic())