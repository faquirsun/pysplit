#!/usr/bin/env python3

"""
The Source class - generate an instance of Source and load stream data and
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

class Source(object):
	"""
	Base Source class
	"""

	def __init__(self, file_path, source_info, station_info):

		self.file_path = file_path

		# Parse out source information
		self.otime = source_info.otime.values[0]
		self.evlat = source_info.evlat.values[0]
		self.evlon = source_info.evlon.values[0]
		self.evdep = source_info.evdep.values[0]
		self.evmag = source_info.evmag.values[0]

		# Parse out station information
		self.slat = station_info.lat.values[0]
		self.slon = station_info.lon.values[0]
		self.selv = station_info.dep.values[0]

		self.sourceid = source_info["sourceid"]

		# Read in source files
		self.stream = read(file_path)

		tmp_stream = self.stream.copy()

		self.component_1 = tmp_stream.select(channel="*Z")[0]
		self.component_2 = tmp_stream.select(channel="*N")[0]
		self.component_3 = tmp_stream.select(channel="*E")[0]

		# Detrend the data
		self.detrendData()

		# Create dictionary to hold window information
		self.window = {'window_beg': None, 'window_end': None}

		# Create dictionary to hold picks
		self.picks = {}

	def loadPicks(self):
		pass

	def loadWindows(self):
		pass

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

	def filterObspy(self, filt_type, minfreq=None, maxfreq=None, n_poles=2, zero_phase=True):
		# Make a copy of the data (so don't have to read a new one if re-filtering)
		tmp_stream = self.stream.copy()

		# Detrend the data
		tmp_stream = (tmp_stream.detrend("linear")).detrend("demean")		

		# Add a cosine taper from the 5th - 95th percentile
		tmp_stream = tmp_stream.taper(max_percentage=0.05)

		# Filter the data
		if filt_type == "Lowpass":
			tmp_stream = tmp_stream.filter(type="lowpass", freq=minfreq, corners=n_poles, zerophase=zero_phase)
		if filt_type == "Bandpass":
			tmp_stream = tmp_stream.filter(type="bandpass", freqmin=minfreq, freqmax=maxfreq, corners=n_poles, zerophase=zero_phase)
		if filt_type == "Highpass":
			tmp_stream = tmp_stream.filter(type="highpass", freq=maxfreq, corners=n_poles, zerophase=zero_phase)

		# Overwrite the 3 components
		self.component_1.data = tmp_stream.select(channel="*Z")[0].data
		self.component_2.data = tmp_stream.select(channel="*N")[0].data
		self.component_3.data = tmp_stream.select(channel="*E")[0].data

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
	
	def _addStat(self, stat, value, pick_type=None):
		if stat == "window_beg" or stat == "window_end":
			self.window[stat] = value

		elif stat == "pick":
			pick = "{}_manual".format(pick_type)
			self.picks[pick] = {"rtime": value,
								"error": self.delta,
								"polarity": "-"} 

		elif stat == "polarity":
			pick = "{}_manual".format(pick_type)
			self.picks[pick]["polarity"] = value

	def saveSource(self, pick_path):
		# Save the traces to mSEED
		file_path_stem = self.file_path[:-1]

		head, tail = os.path.split(file_path_stem)

		# Create a pick file
		for phase, pick in self.picks.items():
			pick_time = (self.starttime + pick["rtime"]).isoformat().split("T")
			ymd = pick_time[0]
			hmsms = pick_time[1].split(":")
			hm = "{}:{}".format(hmsms[0], hmsms[1])
			sms = hmsms[2]
			filename = "{}/{}/{}{}".format(pick_path, self.station.upper(), tail, phase)

			with open(filename, "w+") as f:
				f.write("{} {} {} {} {}".format(ymd, hm, sms, pick["error"], pick["polarity"]))

		for key, value in self.window.items():
			if key == "window_beg":
				extension = "wb"
			elif key == "window_end":
				extension = "we"
			filename = "{}/{}/{}{}".format(pick_path, self.station.upper(), tail, extension)

			with open(filename, "w+") as f:
				f.write("{}".format(value))

	@property
	def distance(self):
		return locations2degrees(self.evlat, self.evlon, self.slat, self.slon)

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