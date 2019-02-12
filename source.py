#!/usr/bin/env python3

"""
Author: Hemmelig
"""

from obspy import UTCDateTime

class Source(object):
	"""
	Base Source class
	"""

	def __init__(self, source_info):

		# Parse out source information
		self._otime     = UTCDateTime(source_info.otime.values[0])
		self._latitude  = source_info.evlat.values[0]
		self._longitude = source_info.evlon.values[0]
		self._depth     = source_info.evdep.values[0]
		self._magnitude = source_info.evmag.values[0]
		self._sourceid  = source_info.sourceid.values[0]

	"""
	Properties from mSEED header:
	            otime : UTCDateTime
	          			Date and time of source origin in UTC
	         latitude : float
						Latitude coordinate of source
	        longitude : float
						Longitude coordinate of source
	          	depth : float
						Depth coordinate of source
	        magnitude : str, optional
						Magnitude of source (Defaults to "-")
	         sourceid : int
	          		    Source ID
	"""

	@property
	def otime(self):
		return self._otime
	
	@property
	def latitude(self):
		return self._latitude
	
	@property
	def longitude(self):
		return self._longitude
	
	@property
	def depth(self):
		return self._depth
	
	@property
	def magnitude(self):
		return self._magnitude
	
	@property
	def sourceid(self):
		return self._sourceid
	