#!/usr/bin/env python3

"""
The Event class - generate an instance of Event and load stream data and
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

class Event(object):
	"""
	Base Event class
	"""

	def __init__(self, file_path, event_info):

		self.file_path = file_path

		# Parse out event information
		self.otime = event_info["otime"]
		self.evlat = event_info["evlat"]
		self.evlon = event_info["evlon"]
		self.evdep = event_info["evdep"]
		self.evmag = event_info["evmag"]

		self.eventid = event_info["sourceid"]

		# Read in event files
		self.stream = read(file_path)

		tmp_stream = self.stream.copy()

		self.Z_comp = tmp_stream.select(channel="*Z")[0]
		self.N_comp = tmp_stream.select(channel="*N")[0]
		self.E_comp = tmp_stream.select(channel="*E")[0]

		# Detrend the data
		self.detrend_data()

		# Create dictionary to hold window information
		self.window = {'window_beg': None, 'window_end': None}

		# Create dictionary to hold picks
		self.picks = {}

	def load_picks(self):
		pass

	def load_windows(self):
		pass

	def unsplit(self, phi, dt):
		# Make copy of the data
		tmp_stream = self.stream.copy()

		# Need to detrend data
		tmp_stream = (tmp_stream.detrend("linear")).detrend("demean")

		# Filter data
		tmp_stream = tmp_stream.filter("bandpass", freqmin=self.filter_min, freqmax=self.filter_max)

		comp_vec = np.array([tmp_stream.select(channel="*N")[0], tmp_stream.select(channel="*E")[0]])

		rot_M = self._make_rot_matrix(phi)

		rot_comp = rot_M @ comp_vec

		lag_comp = [rot_comp[0], np.roll(rot_comp[1], int(dt / self.delta))]

		rot_M = self._make_rot_matrix(-1 * phi)

		unsplit_comp = rot_M @ lag_comp

		return unsplit_comp

	def _make_rot_matrix(self, phi):
		return np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])

	def filter_sac(self, tr, minfreq, maxfreq, n_poles=2, zero_phase=True):

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

	def remove_filter(self):
		self.detrend_data()

	def detrend_data(self):
		# Make a copy of the data
		tmp_stream = self.stream.copy()

		# Detrend the data
		tmp_stream = (tmp_stream.detrend("linear")).detrend("demean")

		# Overwrite the 3 components
		self.Z_comp.data = tmp_stream.select(channel="*Z")[0].data
		self.N_comp.data = tmp_stream.select(channel="*N")[0].data
		self.E_comp.data = tmp_stream.select(channel="*E")[0].data

	def filter_obspy(self, filt_type, minfreq=None, maxfreq=None, n_poles=2, zero_phase=True):
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
		self.Z_comp.data = tmp_stream.select(channel="*Z")[0].data
		self.N_comp.data = tmp_stream.select(channel="*N")[0].data
		self.E_comp.data = tmp_stream.select(channel="*E")[0].data

	def plot_traces(self, Z_ax, N_ax, E_ax, lims):

		if lims == None:
			# Calculate x limits and set them
			xlims = (self.Z_comp.times(type="relative")[0], self.Z_comp.times(type="relative")[-1])
			Z_ax.set_xlim(xlims)
			N_ax.set_xlim(xlims)
			E_ax.set_xlim(xlims)

			# Set y limits to +/- 1
			Z_ax.set_ylim((-1, 1))
			N_ax.set_ylim((-1, 1))
			E_ax.set_ylim((-1, 1))

			# Need to find absolute max of all three traces then normalise them all wrt to that
			z_max = max(abs(self.Z_comp.data))
			n_max = max(abs(self.N_comp.data))
			e_max = max(abs(self.E_comp.data))

			self.norm_factor = max(z_max, n_max, e_max)

		elif lims != None:
			# Set x limits
			xlims = lims[0]
			Z_ax.set_xlim(xlims)
			N_ax.set_xlim(xlims)
			E_ax.set_xlim(xlims)

			# Set y limits
			ylims = lims[1]
			Z_ax.set_ylim(ylims)
			N_ax.set_ylim(ylims)
			E_ax.set_ylim(ylims)

		# Plot the traces
		Z_ax.plot(self.Z_comp.times(type="relative"), self.Z_comp.data / self.norm_factor)
		N_ax.plot(self.N_comp.times(type="relative"), self.N_comp.data / self.norm_factor)
		E_ax.plot(self.E_comp.times(type="relative"), self.E_comp.data / self.norm_factor)

	def _add_stat(self, stat, value, pick_type=None):
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

	def save_event(self, pick_path):
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
		return self.Z_comp.stats.network

	@property
	def station(self):
		return self.Z_comp.stats.station

	@property
	def location(self):
		return self.Z_comp.stats.location
	
	@property
	def channel(self):
		return self.Z_comp.stats.channel
	
	@property
	def starttime(self):
		return self.Z_comp.stats.starttime

	@property
	def endtime(self):
		return self.Z_comp.stats.endtime
	
	@property
	def sampling_rate(self):
		return self.Z_comp.stats.sampling_rate
	
	@property
	def delta(self):
		return self.Z_comp.stats.delta

	@property
	def npts(self):
		return self.Z_comp.stats.npts

	@property
	def calib(self):
		return self.Z_comp.stats.calib