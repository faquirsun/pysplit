#!/usr/bin/env python3

"""
This module requires that 'obspy', 'numpy', 'cartopy' and 'matplotlib' are 
installed in your Python environment.

TO-DO
-----
Come up with a suitable sourceid scheme and implement in _generateCatalogue()

Author: Hemmelig
"""

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
from matplotlib.patches import Circle
import pathlib
import shutil
import metainfo as psm

def parseCatalogueMetafile(meta_dir):
	"""
	Parses the relevant catalogue metafile into a dictionary containing
	the catalogue parameters

	Parameters
	----------
	meta_dir : pathlib.Path object
		Path to metafiles directory

	Returns
	-------

	"""

	params = []

	metafile = meta_dir / "catalogue_metafile.txt"
	with metafile.open(mode='r') as f:
		for line in f:
			if "?" in line:
				line = line.rstrip().split(" ? ")
				params.append(line[1])

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

	return {**catalogue_parameters, **d_specific}


class Catalogue(ABC):
	"""
	Catalogue class

	Abstract Base Class for LocalCatalogue and TeleseismicCatalogue. A catalogue
	contains information pertaining to 

	Attributes
	----------
	src_df : 

	src_cols :

	arr_df :
	
	arr_cols :

	network :



	
	Methods
	-------
	Load catalogue (load a pre-existing catalogue)
	Teleseismic catalogue (create a new teleseismic catalogue - based on functionality from obspy)
	Local catalogue (create a new local catalogue - can be generated from hyp file)
	"""

	src_cols = ['otime', 'lat', 'lon', 'dep', 'mag', 'sourceid']
	arr_cols = ['sourceid', 'receiverid', 'traveltime', 'waveform?']

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
		self.loadSources()
		self.arr_file = self.cat_dir / "metafiles" / "arrivals.txt"

		if new:
			self._generateCatalogueMetafile()
			self._generateCatalogue()

	def getWaveforms(self):
		"""
		Retrieve waveform data from an archive.

		"""

		for i, arrival in self.arr_df.iterrows():
			if arrival['waveform?']:
				continue

			rec = self.network.lookupReceiver(arrival.receiverid)
			src = self.lookupSource(arrival.sourceid)

			self.getWaveform(src, rec, arrival)

			self.arr_df.loc[i, 'waveform?'] = True

			if (i % 10) == 0:
				self.arr_df.to_csv(self.arr_file, index=False)

		self.arr_df.to_csv(self.arr_file, index=False)

	def getWaveform(self, source, receiver, arrival=None):
		"""
		Retrieves waveform data for a single source-receiver pair from an 
		archive.

		Parameters
		----------
		source : Source object

		receiver : Receiver object

		arrival : DataFrame, optional
			DataFrame containing arrival information for the provided source 
			and receiver

		"""

		if not arrival:
			arrival = self.lookupArrival(source.sourceid, receiver.station)

		if arrival["waveform?"].values[0]:
			return True

		data_dir = self.cat_dir / "data" / receiver.station
		data_dir.mkdir(exist_ok=True, parents=True)
		pick_dir = self.cat_dir / "picks" / receiver.station
		pick_dir.mkdir(exist_ok=True, parents=True)
		tmp_dir  = self.cat_dir / "tmp"

		ttime_dict = arrival.traveltime.values[0]

		for phase, ttimes in ttime_dict.items():
			for traveltime in ttimes:
				wbeg = source.otime + float(traveltime) - self.window
				wend = source.otime + float(traveltime) + self.window

				evyear = str(wbeg.year)
				evjday = str(wbeg.julday).zfill(3)
				nxjday = str((wbeg + 86400).julday).zfill(3)

				evfiles = list(self.archive_path.glob(self.archive_format.format(year=evyear, 
																				 jday=evjday, 
																				 receiver=receiver.station, 
																				 comp="*")))
				nxfiles = list(self.archive_path.glob(self.archive_format.format(year=evyear, 
																				 jday=nxjday, 
																				 receiver=receiver.station, 
																				 comp="*")))

				print(evfiles, nxfiles)

				if not evfiles or nxfiles:
					self.arr_df.drop(index=arrival.index.item(), inplace=True)	
					self.arr_df.to_csv(self.arr_file, index=False)

					return False

				for evfile in evfiles:
					shutil.copy(str(evfile), str(tmp_dir))

				if (wbeg.julday != wend.julday):
					for nxfile in nxfiles:
						shutil.copy(str(nxfile), str(tmp_dir))

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
					print("Window required not within time frame of data.")
					for file in files:
						file.unlink()
					continue

				for tr in st:
					name = '{0}/data/{1}/source.{2}.{1}.{3}'.format(self.catalogue_path, 
																	receiver.station, 
															   		arrival.sourceid, 
															   		comp.lower())
					tr.write(name, format="MSEED")

				for file in files:
					file.unlink()

		self.arr_df.at[arrival.index.item(), "waveform?"] = True
		self.arr_df.to_csv(self.arr_file, index=False)

		return True

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
		"""
		Check to see if the waveforms for all of the arrivals in the arrivals
		DataFrame have been downloaded.

		Returns
		-------

		"""

		return self.arr_df["waveform?"].all()

	def lookupSource(self, sourceid):
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

		if type(sourceid) != str:
			print("You must identify the source with a string.")
			return
		try:
			src = self.src_df[self.src_df["sourceid"] == int(sourceid)]
			return psm.Source(src)
		except IndexError:
			print("There is no source with id \"{}\" in this catalogue".format(sourceid))

	def lookupArrival(self, sourceid, receiverid):
		"""
		Queries the arrivals DataFrame for 'sourceid' and 'receiverid'

		Parameters
		----------
		sourceid : str
			Source name/unique identifier
		receiverid : str
			Receiver name/unique identifier

		Returns
		-------
		DataFrame
		
		"""

		return self.arr_df.loc[self.arr_df[["sourceid", "receiverid"]].isin([sourceid, receiverid]).all(1)]

	@property
	def window(self):
		"""Get and set the time window padding either side of the arrival of
		the signal at the receiver. Used when slicing mSEED files. (units: s)"""
		return self._window

	@window.setter
	def window(self, value):
		self._window = value

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


class LocalCatalogue(Catalogue):
	"""
	LocalCatalogue class

	Inherits from Catalogue class, providing functionality specific to 
	catalogues of local sources

	Attributes
	----------
	window : float
		Time window 

	Methods
	-------

	Properties
	----------

	"""

	window = 60.

	def _generateCatalogueMetafile(self):
		"""
		Generates a metafile from a template.

		"""

		filein = open('templates/local_metafile_template.txt')
		src = Template(filein.read())

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

		# Safe substitute writes "" for missing template variables
		output = src.safe_substitute(d_cat)

		meta_path = self.cat_dir / "metafiles" / "catalogue_metafile.txt"
		with meta_path.open(mode='w') as f:
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

		src_path = pathlib.Path(self.local_input)

		if not src_path.exists():
			print("File doesn't exist.")
			print("Please provide a valid path.")
			return

		if self.data_source == "SeisLoc":
			sources = list(src_path.glob("outputs/*.event"))

			# Preallocate memory
			self.src_df = pd.DataFrame(index=np.arange(0, len(sources)), columns=self.src_cols)

			# Needed to track variable lengths of arrival files
			idx = 0
			for i, source in enumerate(sources):
				with source.open(mode='r') as f:
					for j, line in enumerate(f):
						if j == 0:
							continue
						else:
							src_info = line.rstrip().split(",")
							break

				otime = UTCDateTime(src_info[0].replace(" ", "T"))
				lon   = src_info[2]
				lat   = src_info[3]
				dep   = src_info[4]
				mag   = ""
				sourceid = "" ###

				self.src_df.loc[i] = [otime, lat, lon, dep, mag, i]

				arrival = path / "outputs" / "{}.stn".format(source[:-6])
				with arrival.open(mode='r') as f:
					for j, line in enumerate(f):
						if j == 0:
							continue

						line.rstrip().split(",")
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

						self.arr_df.loc[idx] = [i, receiver, ttime, False]

						idx += 1

		elif self.data_source == "PySplit":
			out_dir = self.cat_dir / "metafiles" / src_path.parts[-1]
			shutil.copy(str(src_path), str(out_dir))

			self.src_df = pd.read_csv(input_file)
			self.src_df["otime"] = self.src_df["otime"].apply(UTCDateTime)
			self.src_df["sourceid"] = self.src_df.index

		elif self.data_source == ".hyp file":
			out_dir = self.cat_dir / "metafiles" / src_path.parts[-1]
			shutil.copy(str(src_path), str(out_dir))

			if src_path.suffix == ".hyp":
				lines = []
				with src_path.open(mode='r') as f:
					for line in f:
						if "GEOGRAPHIC" in line:
							lines.append(line)
				lines = [x.rstrip().split() for x in lines]

				# Preallocate memory
				self.src_df = pd.DataFrame(index=np.arange(0, len(lines)), columns=self.src_cols)

				for i, line in enumerate(lines):
					evyear, evmonth, evday  = int(line[2]), int(line[3]), int(line[4])
					evhour, evmin, evsecflt = int(line[5]), int(line[6]), float(line[7])
					evsec  = int(evsecflt)
					evmsec = int(1e6 * (evsecflt - evsec))
					otime  = UTCDateTime(evyear, evmonth, evday, evhour, evmin, evsec, evmsec)

					lat = float(line[9])
					lon = float(line[11])
					dep = float(line[13])
					mag = "?" #### ADD PROPER HANDLING FOR HYP FILES THAT HAVE MAGNITUDES IN
					sourceid = ""

					self.src_df.loc[i] = [otime, lat, lon, dep, mag, i]

			else:
				print("Please provide a .hyp file.")

		elif self.data_source == "IRIS":
			# TO-DO: Implement support for IRIS inputs
			pass

		self.src_df.to_csv(self.src_file, index=False)

	def getArrivals(self, phases=None):
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

		if self.data_source == "PySplit":
			srcs_df = pd.read_csv(self.local_input)
			
			idx = 0
			self.arr_df = pd.DataFrame(index=np.arange(0, len(self.network.receivers) * len(srcs_df)), columns=self.arr_cols)
			for i, source in sources_df.iterrows():
				for j, receiver in self.network.receivers.iterrows():
					# Want to check if the receiver was available at the time of the source
					otime = UTCDateTime(source["otime"])
					stdp  = UTCDateTime(receiver["deployment"])
					etdp  = UTCDateTime(receiver["retrieval"])
					if (otime >= stdp) & (otime <= etdp):
						self.arr_df.loc[idx] = [idx, receiver['receiverid'], "{'{}': [0.0]}", False]
					else:
						self.arr_df.loc[idx] = ["-", "-", "-", "-"]
					idx += 1

		elif self.data_source == ".hyp file":
			srcs = pathlib.Path(self.local_input)
			lines = []
			S_lines = []
			with srcs.open(mode='r') as f:
				for line in f:
					if ">" in line:
						if " P " in line:
							continue
						else:
							lines.append(line)

					if (">" in line) and (" S " in line):
						S_lines.append(line)

			# Preallocate memory
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
					receiver = line[0]
					phase = "S_pred"
					self.arr_df.loc[idx] = [sourceid, receiver, str({phase: [line[15]]}), False]

		# Remove all sources that didn't have complete information
		self.arr_df.drop(self.arr_df[self.arr_df.sourceid == "-"].index, inplace=True)
		self.arr_df.reset_index(drop=True, inplace=True)
		self.arr_df.to_csv(self.arr_file, index=False)

		available_receivers = self.network.filterReceivers(self.arr_df)
		available_receivers.to_csv(self.rec_file, index=False)

	def removeArrival(self, sourceid, receiverid):
		self.arr_df.drop(self.arr_df[(self.arr_df.sourceid == sourceid) \
								   & (self.arr_df.receiverid == receiverid)].index, \
								      inplace=True)
		self.src_df.reset_index(drop=True, inplace=True)
		self.src_df.to_csv(self.src_file, index=False)

	def plotGeographic(self, map_widget, lims=None, receivers=False):
		"""
		Creates a geographic representation of the catalogue.

		Parameters
		----------
		map_widget : Axes object, optional
			The axes object on which to plot the map. One is generated
			if not provided
		lims : list of floats, optional
			Contains the geograhical boundaries of the map. If not provided,
			suitable limits are calculated from the source locations

		"""

		lons = self.src_df.lon.values
		lats = self.src_df.lat.values

		if lims == None: 
			londiff = (max(lons) - min(lons)) * 0.1
			latdiff = (max(lats) - min(lats)) * 0.1
			self.lon0 = min(lons) - londiff
			self.lon1 = max(lons) + londiff
			self.lat0 = min(lats) - latdiff
			self.lat1 = max(lats) + latdiff
		else:
			self.lon0 = lims["lon0"]
			self.lon1 = lims["lon1"]
			self.lat0 = lims["lat0"]
			self.lat1 = lims["lat1"]

		proj = ccrs.PlateCarree()
		self.map = map_widget.canvas.mapPlot(proj)
		self.map.set_xlabel("Longitude, degrees", fontsize=10)
		self.map.set_ylabel("Latitude, degrees", fontsize=10)
		self.map.set_extent([self.lon0, self.lon1, self.lat0, self.lat1], proj)
		self.map.set_aspect('auto')

		self._plotSources()
		if receivers:
			self._plotReceivers()		

	def _plotSources(self):
		lons = self.src_df.lon.values
		lats = self.src_df.lat.values
		sids = self.src_df.sourceid.values

		tolerance = 10
		for i in range(len(lons)):
			self.map.scatter(lons[i], lats[i], 12, marker='o', color='k', picker=tolerance, zorder=10, label="SOURCE: {}".format(sids[i]))

	def _plotReceivers(self):
		lons = self.network.receivers.lon.values
		lats = self.network.receivers.lat.values
		rids = self.network.receivers.name.values
		tolerance = 10
		for i in range(len(lons)):
			self.map.scatter(lons[i], lats[i], 75, marker="v", color="green", picker=tolerance, zorder=15, label="REC: {}".format(rids[i]))


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
		"""
		Generates a metafile from a template.

		"""

		filein = open('templates/teleseismic_metafile_template.txt')
		src = Template(filein.read())

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

		# Safe substitute writes "" for missing template variables
		output = src.safe_substitute(d_cat)

		meta_path = self.cat_dir / "metafiles" / "catalogue_metafile.txt"
		with meta_path.open(mode='w') as f:
			f.write(output)

	def _generateCatalogue(self):
		"""
		Generates a catalogue from the IRIS webclient API.

		"""
		client = Client("IRIS")

		st = UTCDateTime(self.start_date)
		et = UTCDateTime(self.end_date)

		sources = client.get_events(starttime=st, 
									endtime=et, 
									minmagnitude=self.minmag,
									latitude=self.clat, 
									longitude=self.clon, 
									minradius=self.minrad, 
									maxradius=self.maxrad)

		# Preallocate memory
		self.src_df = pd.DataFrame(index=np.arange(0, sources.count()), columns=self.src_cols)

		for i in range(sources.count()):
			try:
				source = sources[i]
				otime = source.preferred_origin().get('time')
				lat = source.preferred_origin().get('latitude')
				lon = source.preferred_origin().get('longitude')
				dep = source.preferred_origin().get('depth') / 1000.0
				mag = source.preferred_magnitude().get('mag')
			except TypeError:
				otime = lat = lon = dep = mag = "-"

			# sourceid = ### Come up with file naming format
			sourceid = i

			self.src_df.loc[i] = [otime, lat, lon, dep, mag, sourceid]

		self.src_df.drop(self.src_df[self.src_df.otime == "-"].index, inplace=True)
		self.src_df.reset_index(drop=True, inplace=True)
		self.src_df.to_csv(self.src_file, index=False)

	def getArrivals(self, phases=["SKS"]):
		"""
		Generates predicted arrival times using a whole earth velocity model

		Parameters
		----------
		phases : list of strings, optional
			List containing the seismic phases for which to calculate 
			predicted traveltimes. Defaults to "SKS"

		"""
		model = TauPyModel(model="ak135")

		rows_list = []
		for i, receiver in self.network.receivers.iterrows():
			rec = psm.Receiver(receiver)
			tmp_df = self.src_df[self.src_df['otime'].between(rec.deployment, rec.retrieval)]

			for j, source in tmp_df.iterrows():
				src = psm.Source(source)
				tmp_dist = locations2degrees(src.latitude, src.longitude, rec.latitude, rec.longitude)
				phase_arr = model.get_travel_times(source_depth_in_km=src.depth,
												   distance_in_degree=tmp_dist,
												   phase_list=phases,
												   receiver_depth_in_km=rec.elevation / 1000)

				if not phase_arr:
					continue
				else:
					phase_dict = {}
					for k in range(len(phase_arr)):
						key  = phase_arr[k].purist_name
						time = phase_arr[k].time
						
						if key in phase_dict:
							phase_dict[key].append(phase_arr[k].time)
						else:
							phase_dict[key] = []
							phase_dict[key].append(time)

							dict1 = {"sourceid": src.sourceid,
									 "receiverid": rec.station,
									 "traveltime": phase_dict,
									 "waveform?": False}
							print(dict1)
							rows_list.append(dict1)

		self.arr_df = pd.DataFrame(rows_list, columns=self.arr_cols)
		self.arr_df.to_csv(self.arr_file, index=False)

		available_receivers = self.network.filterReceivers(self.arr_df)
		available_receivers.to_csv(self.rec_file, index=False)

	def plotGeographic(self, map_widget=None):
		"""
		Creates a geographic representation of the catalogue.

		Parameters
		----------
		map_widget : Axes object, optional
			The axes object on which to plot the map. One is generated
			if not provided

		"""

		network_centre = self.network.centre
		lons = self.src_df.lon.values
		lats = self.src_df.lat.values
		sids = self.src_df.sourceid.values

		proj = ccrs.AzimuthalEquidistant(central_longitude=network_centre[0], central_latitude=network_centre[1])

		if not map_widget:
			return
		else:
			ax = map_widget.canvas.mapPlot(proj)
		
		minr = compute_radius(proj, self.minrad, network_centre)
		ax.add_patch(Circle(xy=network_centre, radius=minr, edgecolor="red",
							linewidth=4.0, fill=False, alpha=0.3, transform=proj,
						    zorder=15))
		if float(self.maxrad) >= 150.:
			pass
		else:
			maxr = compute_radius(proj, self.maxrad, network_centre)
			ax.add_patch(Circle(xy=network_centre, radius=maxr, edgecolor="red",
								linewidth=4.0, fill=False, alpha=0.3, transform=proj, 
								zorder=15))

		tolerance = 10
		for i in range(len(lons)):
			ax.scatter(lons[i], lats[i], 18, marker='o', color='k', picker=tolerance, zorder=20, label="SOURCE: {}".format(sids[i]), transform=ccrs.Geodetic())

def compute_radius(ortho, radius_degrees, netcen):
    phi1 = netcen[1] + float(radius_degrees) if netcen[1] <= 0 else netcen[1] - float(radius_degrees)
    _, y1 = ortho.transform_point(netcen[0], phi1, ccrs.PlateCarree())
    return abs(y1)