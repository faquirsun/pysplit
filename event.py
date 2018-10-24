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

	def __init__(self, file_path):

		# Read in event files
		self.stream = read(file_path)

		# Select the 3 components
		self.Z_comp = self.stream.select(channel="*Z")[0]
		self.N_comp = self.stream.select(channel="*N")[0]
		self.E_comp = self.stream.select(channel="*E")[0]

		# Add window dictionary to each component
		self.Z_comp.stats.window = {'window_beg': None, 'window_end': None}
		self.N_comp.stats.window = {'window_beg': None, 'window_end': None}
		self.E_comp.stats.window = {'window_beg': None, 'window_end': None}

		# Add empty picks dictionary to each component
		self.Z_comp.stats.picks = {}
		self.N_comp.stats.picks = {}
		self.E_comp.stats.picks = {}

		self.delta = self.Z_comp.stats.delta

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
		self.Z_comp = self.stream.select(channel="*Z")[0]
		self.N_comp = self.stream.select(channel="*N")[0]
		self.E_comp = self.stream.select(channel="*E")[0]

	def filter_obspy(self, filt_type, minfreq, maxfreq, n_poles=2, zero_phase=True):
		# Make a copy of the data (so don't have to read a new one if re-filtering)
		tmp_stream = self.stream.copy()

		# Detrend the data
		tmp_stream = (tmp_stream.detrend("linear")).detrend("demean")

		# Add a cosine taper from the 5th - 95th percentile
		tmp_stream = tmp_stream.taper(max_percentage=0.05)

		# Filter the data
		tmp_stream.filter(type=filt_type, freqmin=minfreq, freqmax=maxfreq, corners=n_poles, zerophase=zero_phase)

		# Overwrite the 3 components
		self.Z_comp = tmp_stream.select(channel="*Z")[0]
		self.N_comp = tmp_stream.select(channel="*N")[0]
		self.E_comp = tmp_stream.select(channel="*E")[0]

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
			self.Z_comp.stats.window[stat] = value
			self.N_comp.stats.window[stat] = value
			self.E_comp.stats.window[stat] = value

		elif stat == "pick":
			pick = "{}_manual".format(pick_type)
			self.Z_comp.stats.picks[pick] = value
			self.N_comp.stats.picks[pick] = value
			self.E_comp.stats.picks[pick] = value