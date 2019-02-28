#!/usr/bin/env python3

"""
This module requires that 'obspy' and 'numpy' are installed in your Python
environment.

TO-DO
-----
Consider breaking the 'Catalogue' class into 'Network' and 'Catalogue' and
an object that describes the matching between these similar to 
SourceReceiverPair

Author: Hemmelig
"""

from obspy import read, UTCDateTime
from obspy.geodetics import locations2degrees, gps2dist_azimuth
import numpy as np
import ast
import pathlib
import pandas as pd

class SourceReceiverPair(object):
	"""
	Object SourceReceiver 

	This class is used to represent the waveforms of a particular event recorded
	at a particular station.

	Attributes
	----------
	file_format : str
		formattable string used for naming output files
	source : Source-type object
		a Source object that contains information specific to the source
	receiver : Receiver-type object
		a Receiver object that contains information specific to the receiver
	components : str
		a string that tracks the current component system ("ZNE" or "ZRT")
	pick_path : str
		a formatted string that points to the location of the relevant pick
		files

	Methods
	-------
	filter(code=None, filt)
		Filters the seismic data
	rotate(method)
		Rotates the seismic data
	saveData()
		Saves manual pick data
	addData(data, value, pick_type=None)
		Adds pick information to dictionary

	Properties
	----------
	network : string, optional
		Network code (default is an empty string)
	station : string, optional
		Station code
	location : string, optional
		Location code (default: "")
	starttime : UTCDateTime, optional
		Date and time of the first data sample given in UTC
	endtime : UTCDateTime, optional
		Date and time of last data sample given in UTC
	sampling_rate : float, optional
		Sampling rate (units: Hz; default value: 1.0)
	delta : float, optional
		Sample distance (units: s; default value: 1.0)
	npts : int, optional
		Number of sample points (default value: 0)
	calib : float, optional
		Calibration factor (default value: 1.0) 

	"""

	file_format = "source.{}.{}"

	def __init__(self, catalogue_path, source, receiver):
		"""
		Object initialisation

		Parameters
		----------
		catalogue_path : str
			Path to catalogue containing this source-receiver pair
		source : Source object
			a Source object that contains information specific to the source
		receiver : Receiver object
			a Receiver object that contains information specific to the receiver
		"""

		self.source   = source
		self.receiver = receiver

		# Set source-receiver filename
		srp_fname = self.file_format.format(source.sourceid, receiver.station)

		# Pick files
		pick_dir = pathlib.Path(catalogue_path) / "picks" / receiver.station
		pick_dir.mkdir(exist_ok=True)
		self.pick_path = pick_dir / srp_fname

		# Source files
		data_dir = pathlib.Path(catalogue_path) / "data" / receiver.station
		src_path = data_dir / srp_fname

		self.stream = read(source_file)

		self.components = "ZNE"

		self._distance, self._az, self._baz = gps2dist_azimuth(source.latitude, source.longitude, 
															   receiver.latitude, receiver.longitude)

		self._updateComponents(self.stream)
		self._loadPicks()

	def filter(self, method=None, filt=None):
		"""
		Filters the seismic data and updates the component attributes.

		Parameters
		----------
		method : str, optional
			The method to be used for filtering. Options are "obspy", "sac" or
			None (default to None, which removes the filter)
		filt : dict
			Contains information about the filter
		"""

		tmp_stream = self.stream.copy()

		if code == "sac":
			### TO-DO : Adapt to handle filter in dictionary format
			###			Move to pathlib from os
			pass
			# if filt["zero_phase"]:
			# 	ph = 2
			# else:
			# 	ph = 1
			# tmp_txt = """
			#			sac<< EOF > /dev/null\n
			#  			r f.sac\n
			# 			rmean\n
			#			taper\n
			#	  		bp c {:f} {:f} n {:d} p {:d} \n
			#			write over\n
			#			q\n
			#			EOF\n
			#			"""
			# tmp_dir = pathlib.(os.path.expanduser('~'), 'filter_temp')
			# if not(os.path.exists(temp_dir)):
			# 	os.makedirs(temp_dir)
			# tr.write('{}/f.sac'.format(temp_dir), "SAC")
			# text = temp_text.format(minfreq, maxfreq, n_poles, ph)
			# os.chdir(temp_dir)
			# open('s', 'w+').write(text)
			# os.system('bash s')
			# tr = read('{}/f.sac'.format(temp_dir))[0]
			# return tr
		else:
			(tmp_stream.detrend("linear")).detrend("demean")	
			if code == "obspy":
				tmp_stream.taper(max_percentage=0.05)
				tmp_stream.filter(**filt)

		self._updateComponents(tmp_stream)

	def rotate(self, method=None, phi=None):
		"""Rotates the seismic data and updates the component attributes.

		Parameters
		----------
		method : str
			Specifies the rotation operation to be used
		phi : float
			Angle through which to rotate components (units: degrees)
		"""
		if not method:
			rphi = np.deg2rad(phi)

		elif method == "NE->RT":
			self.components = "ZRT"
		elif method == "RT->NE":
			self.components = "ZNE"

		self.stream.rotate(method, back_azimuth=self.baz)

		self._updateComponents(tmp_stream)

	def saveData(self):
		"""Save pick data.

		"""
		p = self.pick_path.with_ext(".pf")
		with open(p, 'w') as f:
			print(self.picks, file=f)

	def addPick(self, info_type, value, pick_type):
		"""Add pick information to dictionary

		Parameters
		----------
		info_type : str
			Specifies what type of pick information is being added			
		value : float
			Time of pick relative to start time of trace (unit: s)
		pick_type : str
			Specifies what type of pick it is - either a specific phase
			pick or a time window
		"""
		if info_type == "pick":
			self.picks[pick_type] = {"rtime": value,
									 "error": self.delta,
									 "polarity": "-"}

		elif info_type == "polarity":
			self.picks[pick_type][info_type] = value

		elif info_type == "window":
			self.picks[pick_type] = {"rtime": value}

	def _updateComponents(self, stream):
		"""Update the component attributes
		
		Parameters
		----------
		stream : Stream object
			List like object of multiple ObsPy Trace objects
		"""

		if self.components == "ZNE":
			self.component_1 = stream.select(channel="*Z")[0]
			self.component_2 = stream.select(channel="*N")[0]
			self.component_3 = stream.select(channel="*E")[0]
		elif self.components == "ZRT":
			self.component_1 = stream.select(channel="*Z")[0]
			self.component_2 = stream.select(channel="*R")[0]
			self.component_3 = stream.select(channel="*T")[0]

	def _loadPicks(self):
		"""Load pick data.

		"""
		p = self.pick_path.with_ext(".pf")
		try:
			with p.open(pick_file) as f:
				picks = f.readline()
				picks = ast.literal_eval(picks)
		except FileNotFoundError:
			print("Currently no pick file.")
			picks = {}

		self._picks = picks

	# --------------------------

	def correctSplitting(self, phi, dt):
		# Make copy of the data
		tmp_stream = self.stream.copy()

		# Need to detrend data
		tmp_stream = (tmp_stream.detrend("linear")).detrend("demean")

		# Filter data
		tmp_stream = tmp_stream.filterObspy(self.filt)

		# Create vector of N and E components and correct splitting
		components = np.array([tmp_stream.select(channel="*N")[0], 
							   tmp_stream.select(channel="*E")[0]])
		rot_M = self._makeRotationMatrix(phi)
		rotated_components = rot_M @ components
		lagged_components = [rotated_components[0], np.roll(rotated_components[1], int(dt / self.delta))]
		rot_M = self._makeRotationMatrix(-1 * phi)
		corrected_components = rot_M @ lagged_components

		return corrected_components[0], corrected_components[1]

	### TO-DO : Add handling for when no axes are provided (should be default,
	###			as it makes it stand-alone from seispick interface).
	def plotTraces(self, axes=None, lims=None):
		"""Plot component traces on provided axes with specified limits

		Parameters
		----------
		axes : list of Axes objects, optional
			An array containing the axes
		lims : list of float lists

		"""

		if lims == None:
			xlims = (self.component_1.times(type="relative")[0], 
					 self.component_1.times(type="relative")[-1])
			ylims = (-1, 1)

		elif lims != None:
			xlims = lims[0]
			ylims = lims[1]

		axes[0].set_xlim(xlims)
		axes[1].set_xlim(xlims)
		axes[2].set_xlim(xlims)

		axes[0].set_ylim(ylims)
		axes[1].set_ylim(ylims)
		axes[2].set_ylim(ylims)

		axes[0].plot(self.component_1.times(type="relative"), self.component_1.data / self.norm_factor, linewidth=0.5)
		axes[1].plot(self.component_2.times(type="relative"), self.component_2.data / self.norm_factor, linewidth=0.5)
		axes[2].plot(self.component_3.times(type="relative"), self.component_3.data / self.norm_factor, linewidth=0.5)

		axes[0].set_ylabel("Z")
		if self.components == "ZNE":
			axes[1].set_ylabel("N")
			axes[2].set_ylabel("E")

		elif self.components == "ZRT":
			axes[1].set_ylabel("R")
			axes[2].set_ylabel("T")

	def plotSpectrogram(self, ax=None, log=True):
		"""Plot spectrogram of the stream data

		Parameters
		----------
		ax : Axes object, optional
			Axes object on which the spectrogram is plotted
		log : Boolean
			Boolean variable that specifies whether to use a log scale
			on the y-axis of the spectrogram (default True)
		"""

		self.stream.spectrogram(log=log, axes=ax)
		ax.set_ylabel("log(Power)")

	def plotHodogram(self, ax, phi=None, dt=None):
		### TO DO: Grab the windows (if they exist)

		if not phi and not dt:
			x = self.component_3.data
			y = self.component_2.data
		else:
			x, y = self.correctSplitting(phi, dt)

		# Normalise the particle motion plots
		norm = max(max(abs(x)), max(abs(y)))

		x = [z / norm for z in x]
		y = [z / norm for z in y]

		# Plot to axes
		ax.plot(x, y)

	def SNR(self, method=None, phase=None):
		# Check if phase and method have been provided
		if not phase:
			return
		if not method:
			return

	def sourcePolarisation(self):
		pass

	@property
	def picks(self):
		"""Get and set values in the dictionary that tracks picks
		related to this source-receiver pair"""
		return self._picks

	@picks.setter
	def picks(self, value):
		self._picks = value

	@property
	def components(self):
		"""Get and set the component coordinate system"""
		return self._components

	@components.setter
	def components(self, value):
		self._components = value	

	@property
	def distance(self):
		"""Get the distance in epicentral angle between source and receiver"""
		return self._distance

	@property
	def baz(self):
		"""Get the back-azimuth between source and receiver"""
		return self._baz

	@property
	def az(self):
		"""Get the azimuth between source and receiver"""
		return self._az

	@property
	def norm_factor(self):
		"""Get the normalising factor for the traces"""
		z_max = max(abs(self.component_1.data))
		n_max = max(abs(self.component_2.data))
		e_max = max(abs(self.component_3.data))

		return max(z_max, n_max, e_max)
	
	@property
	def network(self):
		"""Get the network code"""
		return self.receiver.network

	@property
	def station(self):
		"""Get the station code"""
		return self.receiver.station

	@property
	def location(self):
		"""Get the geographical coordinates of the source and receiver"""
		return self.source.lon, self.source.lat, self.source.dep, \
			   self.receiver.lon, self.receiver.lat, self.receiver.elv
	
	@property
	def starttime(self):
		"""Get the date and time of the first data sample"""
		return self.component_1.stats.starttime

	@property
	def endtime(self):
		"""Get the date and time of the last data sample"""
		return self.component_1.stats.endtime
	
	@property
	def sampling_rate(self):
		"""Get the sampling rate"""
		return self.component_1.stats.sampling_rate
	
	@property
	def delta(self):
		"""Get the sample distance"""
		return self.component_1.stats.delta

	@property
	def npts(self):
		"""Get the number of data"""
		return self.component_1.stats.npts

	@property
	def calib(self):
		"""Get the calibration factor"""
		return self.component_1.stats.calib


class Receiver(object):
	"""
	Receiver class

	Object to capture the information that uniquely describes a receiver,
	including the network/station identifier codes, geographical location and 
	deployment period.

	Properties
	----------
	network : str, optional
		The network code (default is an empty string)
	station : str
		The station code
	latitude : float
		Latitude coordinate of receiver
	longitude : float
		Longitude coordinate of receiver
	elevation : float
		Elevation coordinate of receiver (unit: m)
	deployment : UTCDateTime, optional
		Date and time of receiver deployment given in UTC
		(default is an empty string)
	retrieval : UTCDateTime, optional
		Date and time of receiver retrieval given in UTC
		(default is an empty string)
	"""

	def __init__(self, rec_info):
		"""
		Parameters
		----------
		rec_info : pandas DataFrame
			Contains the information to initialise a receiver object
		"""

		try:
			self._network = rec_info.network.values[0]
		except AttributeError:
			self._network = ""
		self._station   = rec_info.name.values[0]
		self._latitude  = rec_info.lat.values[0]
		self._longitude = rec_info.lon.values[0]
		self._elevation = rec_info.elv.values[0]
		try:
			self._deployment = UTCDateTime(rec_info.deployment.values[0])
			self._retrieval  = UTCDateTime(rec_info.retrieval.values[0])
		except AttributeError:
			self._deployment = ""
			self._retrieval  = ""

	@property
	def network(self):
		"""Get the network code"""
		return self._network
	
	@property
	def station(self):
		"""Get the station code"""
		return self._station.upper()

	@property
	def latitude(self):
		"""Get the receiver latitude"""
		return self._latitude

	@property
	def longitude(self):
		"""Get the receiver longitude"""
		return self._longitude

	@property
	def elevation(self):
		"""Get the receiver elevation"""
		return self._elevation

	@property
	def deployment(self):
		"""Get the receiver deployment date as a UTCDateTime object"""
		return self._deployment

	@property
	def retrieval(self):
		"""Get the receiver retrieval date as a UTCDateTime object"""
		return self._retrieval
		

class Source(object):
	"""Source class

	Object to capture the information that uniquely describes a source,
	including the origin time, hypocenter location and magnitude.

	Properties
	----------
	otime : UTCDateTime
		Date and time of source origin time given in UTC
	latitude : float
		Latitude coordinate of source
	longitude : float
		Longitude coordinate of source
	depth : float
		depth coordinate of source (unit: km)
	magnitude : str, optional
		magnitude of source
	"""

	def __init__(self, src_info):
		"""
		Parameters
		----------
		src_info : pandas DataFrame
			Contains the information to initialise a source object
		"""

		self._otime     = UTCDateTime(src_info.otime.values[0])
		self._latitude  = src_info.lat.values[0]
		self._longitude = src_info.lon.values[0]
		self._depth     = src_info.dep.values[0]
		self._magnitude = src_info.mag.values[0]
		self._sourceid  = src_info.sourceid.values[0]

	@property
	def otime(self):
		"""Get the origin time"""
		return self._otime
	
	@property
	def latitude(self):
		"""Get the source latitude""" 
		return self._latitude
	
	@property
	def longitude(self):
		"""Get the source longitude"""
		return self._longitude
	
	@property
	def depth(self):
		"""Get the source depth"""
		return self._depth
	
	@property
	def magnitude(self):
		"""Get and set the source magnitude. Returns None if no magnitude"""
		return self._magnitude

	@magnitude.setter
	def magnitude(self, value):
		if type(value) == float:
			self._magnitude = value
		else:
			#raise TypeError("Magnitudes must be of type 'float'")
			self._magnitude = None


class Network(object):
	"""
	Network class

	Attributes
	----------
	rec_cols : list
		Contains header strings for receiver DataFrame
	receivers : DataFrame
		Pandas DataFrame containing receiver information

	Methods
	-------
	lookupReceiver(receiverid)
		Looks up and returns a Receiver object for a given string identifier

	Examples
	--------
	>>> import pysplit.metainfo as psm
	>>> nw = psm.Network("~/path/to/file")
	>>> b = nw.lookupReceiver("rec1")
	>>> b.latitude
	65.42

	"""

	rec_cols = ['network', 'name', 'lat', 'lon', 'elv', 'deployment', 'retrieval']

	def __init__(self, rec_file):
		"""
		Class initialisation

		Parameters
		----------
		rec_file : pathlib Path object
			Path pointing to csv file containing receiver information

		"""
		super().__init__()

		self.receivers = pd.DataFrame(columns=self.rec_cols)

		if type(rec_file) == str:
			rec_file = pathlib.Path(rec_file)

		self._loadReceivers(rec_file)

	def lookupReceiver(self, receiverid):
		"""
		Queries the network DataFrame for 'receiverid'.

		Parameters
		----------
		receiverid : str
			Receiver name/unique identifier

		Returns
		-------
		Receiver object
		"""
		if not type(receiverid) == str:
			print("You must identify the receiver with a string.")
			return
		try:
			rec = self.receivers.loc[self.receivers["name"] == receiverid]
			return Receiver(rec)
		except IndexError:
			print("There is no receiver with id \"{}\" in this network".format(receiverid))

	def filterReceivers(self, arrivals):
		available = arrivals.drop_duplicates(subset="receiverid", keep="first").name.values

		return self.receivers.loc[self.receivers['name'].isin(available)]

	def _loadReceivers(self, rec_file):
		if rec_file.is_file():
			self.receivers = pd.read_csv(rec_file)
			try:
				self.receivers['deployment'] = self.receivers['deployment'].apply(UTCDateTime)
				self.receivers['retrieval']  = self.receivers['retrieval'].apply(UTCDateTime)
			except TypeError:
				pass

	@property
	def centre(self):
		"""Get the geographical centre of the network."""
		return [self.rec_df.lon.mean(), self.rec_df.lat.mean()]