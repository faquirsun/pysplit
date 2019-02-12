#!/usr/bin/env python3

"""
The SourceReceiverPair class - generate an instance of SourceReceiverPair and load stream data and
metadata.

Methods
Split: 

Attributes

Author: Hemmelig
"""

import obspy
from obspy import read
import numpy as np
import os
from obspy import UTCDateTime
from obspy.geodetics import locations2degrees
import ast

class SourceReceiverPair(object):
	"""
	Base Source class
	"""

	# --------------------
	# Class initialisation
	# --------------------

	def __init__(self, catalogue_path, source, receiver):

		self.source   = source
		self.receiver = receiver

		file_path = "{0}/data/{1}/source.{2}.{1}.*".format(catalogue_path, self.receiver.station.upper(), self.source.sourceid)
		self.pick_path = "{0}/picks/{1}/source.{2}.{1}".format(catalogue_path, self.receiver.station.upper(), self.source.sourceid)
		if not os.path.exists(self.pick_path):
			os.makedirs(self.pick_path)

		# Read in source files
		self.stream = read(file_path)

		tmp_stream = self.stream.copy()

		self.component_1 = tmp_stream.select(channel="*Z")[0]
		self.component_2 = tmp_stream.select(channel="*N")[0]
		self.component_3 = tmp_stream.select(channel="*E")[0]

		# Detrend the data
		self.detrendData()

		# Lookup window picks
		self._window = self.loadWindows()

		# Lookup phase picks
		self._picks = self.loadPicks()

	# --------------------

	# --------------------------
	# Functions related to picks 
	# --------------------------

	def loadWindows(self):
		# Define path to window file (if exists)
		window_file = "{}.wf".format(self.pick_path)
		# Read in window file
		try:
			with open(window_file) as f:
				window = f.readline()
				window = ast.literal_eval(window)
		except FileNotFoundError:
			print("Currently no window file.")
			window = {}

		return window

	def loadPicks(self):
		# Define path to pick file (if exists)
		pick_file = "{}.pf".format(self.pick_path)
		# Read in pick file
		try:
			with open(pick_file) as f:
				picks = f.readline()
				picks = ast.literal_eval(picks)
		except FileNotFoundError:
			print("Currently no pick file.")
			picks = {}

		return picks

	def saveData(self):
		pick_file = "{}.pf".format(self.pick_path)
		with open(pick_file, 'w') as f:
			print(self.picks, file=f)

		window_file = "{}.wf".format(self.pick_path)
		with open(window_file, 'w') as f:
			print(self.window, file=f)

	def addData(self, data, value, pick_type=None):
		if data == "pick":
			pick = "{}_manual".format(pick_type)
			self._picks[pick] = {"rtime": value,
								"error": self.delta,
								"polarity": "-"} 

		elif data == "polarity":
			pick = "{}_manual".format(pick_type)
			self._picks[pick]["polarity"] = value

		elif data == "wb" or data == "we":
			self._window[data] = value

	# --------------------------

	def unsplit(self, phi, dt):
		# Make copy of the data
		tmp_stream = self.stream.copy()

		# Need to detrend data
		tmp_stream = (tmp_stream.detrend("linear")).detrend("demean")

		# Filter data
		tmp_stream = tmp_stream.filter("bandpass", freqmin=self.filter_min, freqmax=self.filter_max)

		comp_vec = np.array([tmp_stream.select(channel="*N")[0], tmp_stream.select(channel="*E")[0]])

		rot_M = self._makeRotationMatrix(phi)

		rot_comp = rot_M @ comp_vec

		lag_comp = [rot_comp[0], np.roll(rot_comp[1], int(dt / self.delta))]

		rot_M = self._makeRotationMatrix(-1 * phi)

		unsplit_comp = rot_M @ lag_comp

		return unsplit_comp

	def _makeRotationMatrix(self, phi):
		return np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])

	# ----------------
	# Filter functions
	# ----------------

	def filterObspy(self, filt):
		# Make a copy of the data (so don't have to read a new one if re-filtering)
		tmp_stream = self.stream.copy()

		# Detrend the data
		tmp_stream = (tmp_stream.detrend("linear")).detrend("demean")		

		# Add a cosine taper from the 5th - 95th percentile
		tmp_stream = tmp_stream.taper(max_percentage=0.05)

		# Filter the data
		tmp_stream = tmp_stream.filter(**filt)

		# Overwrite the 3 components
		self.component_1.data = tmp_stream.select(channel="*Z")[0].data
		self.component_2.data = tmp_stream.select(channel="*N")[0].data
		self.component_3.data = tmp_stream.select(channel="*E")[0].data

	### TO-DO: Adapt to handle filter in dictionary format
	def filterSAC(self, tr, minfreq, maxfreq, n_poles=2, zero_phase=True):

		if zero_phase:
			p = 2
		else:
			p = 1
	
		temp_text = """sac<< EOF > /dev/null\nr f.sac\nrmean\ntaper\nbp c {:f} {:f} n {:d} p {:d} \nwrite over\nq\nEOF\n"""
		temp_dir = os.path.join(os.path.expanduser('~'), 'filter_temp')
		if not(os.path.exists(temp_dir)):
			os.makedirs(temp_dir)

		tr.write('{}/f.sac'.format(temp_dir), "SAC")
		text = temp_text.format(minfreq, maxfreq, n_poles, p)
		os.chdir(temp_dir)
		open('s', 'w+').write(text)
		os.system('bash s')
		tr = obspy.read('{}/f.sac'.format(temp_dir))[0]
		return tr

	def removeFilter(self):
		# Remove the data by simply calling the detrend data function
		self.detrendData()

	def detrendData(self):
		# Make a copy of the data
		tmp_stream = self.stream.copy()

		# Detrend the data
		tmp_stream = (tmp_stream.detrend("linear")).detrend("demean")

		# Overwrite the 3 components
		self.component_1.data = tmp_stream.select(channel="*Z")[0].data
		self.component_2.data = tmp_stream.select(channel="*N")[0].data
		self.component_3.data = tmp_stream.select(channel="*E")[0].data

	# ----------------

	# ------------------
	# Plotting functions
	# ------------------

	def plotTraces(self, component_1_axis, component_2_axis, component_3_axis, lims):

		if lims == None:
			# Calculate x limits and set them
			xlims = (self.component_1.times(type="relative")[0], self.component_1.times(type="relative")[-1])
			component_1_axis.set_xlim(xlims)
			component_2_axis.set_xlim(xlims)
			component_3_axis.set_xlim(xlims)

			# Set y limits to +/- 1
			component_1_axis.set_ylim((-1, 1))
			component_2_axis.set_ylim((-1, 1))
			component_3_axis.set_ylim((-1, 1))

			# Need to find absolute max of all three traces then normalise them all wrt to that
			z_max = max(abs(self.component_1.data))
			n_max = max(abs(self.component_2.data))
			e_max = max(abs(self.component_3.data))

			self.norm_factor = max(z_max, n_max, e_max)

		elif lims != None:
			# Set x limits
			xlims = lims[0]
			component_1_axis.set_xlim(xlims)
			component_2_axis.set_xlim(xlims)
			component_3_axis.set_xlim(xlims)

			# Set y limits
			ylims = lims[1]
			component_1_axis.set_ylim(ylims)
			component_2_axis.set_ylim(ylims)
			component_3_axis.set_ylim(ylims)

		# Plot the traces
		component_1_axis.plot(self.component_1.times(type="relative"), self.component_1.data / self.norm_factor, linewidth=0.5)
		component_2_axis.plot(self.component_2.times(type="relative"), self.component_2.data / self.norm_factor, linewidth=0.5)
		component_3_axis.plot(self.component_3.times(type="relative"), self.component_3.data / self.norm_factor, linewidth=0.5)

	def plotSpectrogram(self, ax):
		self.stream.spectrogram(log=True, axes=ax)

	# ------------------

	@property
	def picks(self):
		return self._picks

	@property
	def window(self):
		return self._window

	@property
	def distance(self):
		return locations2degrees(self.source.latitude, self.source.longitude, self.receiver.latitude, self.receiver.longitude)

	"""
	Properties from mSEED header:
	          network : string, optional
	          			Network code (default is an empty string)
	          station : string, optional
	          			Station code (default is an empty string)
	         location : string, optional
	         			Location code (default is an empty string)
	          channel : string, optional
	          			Channel code (default is an empty string)
		    starttime : UTCDateTime, optional
		    			Date and time of the first data sample given in UTC
		      endtime : UTCDateTime, optional
		      			Date and time of last data sample given in UTC
		sampling_rate : float, optional
						Sampling rate in hertz (default value is 1.0)
				delta : float, optional
						Sample distance in seconds (default value is 1.0)
				 npts : int, optional
				 		Number of sample points (default value is 0, which implies no data is present)
				calib : float, optional
						Calibration factor (default value is 1.0) 
	"""
	@property
	def network(self):
		return self.component_1.stats.network

	@property
	def station(self):
		return self.component_1.stats.station

	@property
	def location(self):
		return self.component_1.stats.location
	
	@property
	def channel(self):
		return self.component_1.stats.channel
	
	@property
	def starttime(self):
		return self.component_1.stats.starttime

	@property
	def endtime(self):
		return self.component_1.stats.endtime
	
	@property
	def sampling_rate(self):
		return self.component_1.stats.sampling_rate
	
	@property
	def delta(self):
		return self.component_1.stats.delta

	@property
	def npts(self):
		return self.component_1.stats.npts

	@property
	def calib(self):
		return self.component_1.stats.calib