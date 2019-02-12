#!/usr/bin/env python3

"""
Author: Hemmelig
"""

from obspy import UTCDateTime

class Receiver(object):
	"""
	Base Receiver class
	"""

	def __init__(self, receiver_info):

		# Parse out receiver information
		try:
			self._network = receiver_info.network.values[0]
		except AttributeError:
			print("No network in receiver file")
		self._station = receiver_info.receiver_name.values[0]
		self._latitude   = receiver_info.lat.values[0]
		self._longitude  = receiver_info.lon.values[0]
		self._elevation  = receiver_info.dep.values[0]
		self._start_deployment = UTCDateTime(receiver_info.st_dep.values[0])
		self._end_deployment   = UTCDateTime(receiver_info.et_dep.values[0])
		self._receiverid = receiver_info.receiverid.values[0]

	"""
	Receiver properties:
	          network : string
	          			Network code (default is an empty string)
	          station : string
	          			Station code (default is an empty string)
	         latitude : float
	          			Latitude coordinate of receiver
		    longitude : float
		    			Longitude coordinate of receiver
		    elevation : float
		      			Elevation coordinate of receiver
		   start_date : UTCDateTime, optional
		      			Date and time of last data sample given in UTC
			 end_date : UTCDateTime, optional
		      			Date and time of last data sample given in UTC
		   receiverid : int
				 		Receiver ID
	"""

	@property
	def network(self):
		return self._network
	
	@property
	def station(self):
		return self._station

	@property
	def latitude(self):
		return self._latitude

	@property
	def longitude(self):
		return self._longitude

	@property
	def elevation(self):
		return self._elevation

	@property
	def start_deployment(self):
		return self._start_deployment

	@property
	def end_deployment(self):
		return self._end_deployment

	@property
	def receiverid(self):
		return self._receiverid