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

		self.stream = read(file_path)

		self.Z_comp = self.stream.select(channel="*Z")[0]
		self.N_comp = self.stream.select(channel="*N")[0]
		self.E_comp = self.stream.select(channel="*E")[0]

		# Need to find filter parameters
		self.filter_min, self.filter_max = (1, 10)

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

	def filter_with_sac(self, tr, minfreq, maxfreq, n_poles=2, zero_phase=True):

		if zero_phase:
			p = 2
		else:
			p = 1
	
		temp_text = """sac<< EOF > /dev/null\nr f.sac\nrmean\ndif\nmul -1\ntaper\nbp c {:f} {:f} n {:d} p {:d} \nwrite over\nq\nEOF\n"""
		temp_dir = os.path.join(os.path.expanduser('~'), 'filter_temp')
		if not(os.path.exists(temp_dir)):
			os.makedirs(temp_dir)

		tr.write('{}/f.sac'.format(temp_dir), "SAC")
		text = temp_text.format(minfreq, maxfreq, n_poles, p)
		print(text)
		os.chdir(temp_dir)
		open('s', 'w+').write(text)
		os.system('bash s')
		tr = obspy.read('{}/f.sac'.format(temp_dir))[0]
		return tr

	def plot_traces(self, Z_ax, N_ax, E_ax):

		Z_ax.plot(self.Z_comp.times(type="relative"), self.Z_comp.data)
		N_ax.plot(self.N_comp.times(type="relative"), self.N_comp.data)
		E_ax.plot(self.E_comp.times(type="relative"), self.E_comp.data)