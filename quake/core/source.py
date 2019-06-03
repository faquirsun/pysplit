# -*- coding: utf-8 -*-
"""
This module requires that "obspy" is installed in your Python environment.

Author: Conor Bacon
"""

from obspy import UTCDateTime

from quake.core.arrival import Arrival
from quake.core.pick import Pick


class Source(object):
    """
    Source class

    The Source class encapsulates the information that uniquely describes a
    seismic event, including origin time, hypocentral location, magnitude and
    a composition of Arrival objects associated with the Source.

    Attributes
    ----------
    uid : str
        The source code
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
    arrivals: dict of Arrival objects
        Composition of Arrivals associated with this Source

    Methods
    -------
    lookup(receiverid)
        Queries the source with the identifier and returns an Arrival object
    get_iris_arrivals(phases, model, receivers)
        Populates Source with Arrival objects using obpsy TauPyModel
    get_hyp_arrival(receiver, pick)
        Parses a line from a .hyp file into an Arrival object
    load_waveforms()
        Read waveform data into all Arrivals associated with this Source

    """

    def __init__(self, src_info, network=None):
        """
        Class initialisation

        Parameters
        ----------
        src_info : pandas DataFrame
            Contains the information to initialise a source object

        """

        super().__init__()

        self.network = network

        for key, value in src_info.items():
            setattr(self, key, value)

        if isinstance(self.otime, str):
            self.otime = UTCDateTime(self.otime)

        if not hasattr(self, "arrivals"):
            self.arrivals = {}

    def lookup(self, receiverid):
        """
        Queries the source for a given arrival

        Parameters
        ----------
        receiverid : str
            Receiver name/unique identifier

        Returns
        -------
        Arrival object

        """

        if receiverid in self.arrivals:
            arrival = self.arrivals[receiverid]
            return arrival
        else:
            msg = "There is no arrival recorded at '{}' for this source."
            print(msg.format(receiverid))

    def load(self):
        """
        Load the arrivals in self.arrivals into Arrival objects

        """

        arrivals = self.arrivals.copy()

        self.arrivals = {}

        for uid, arrival in arrivals.items():
            arr = Arrival(self, self.network.lookup(uid))
            arr.load_picks(arrival["picks"])
            self + arr

    def get_iris_arrivals(self, phases, model, receivers):
        """
        Get the predicted traveltimes of the specified seismic phases at the
        set of receivers and generate a suite of Arrival objects

        Parameters
        ----------
        phases : list-like
            List of seismic phase codes for which to calculate arrival times
        model : TauPyModel object
            Velocity model within which to calculate arrival times
        receivers : dict
            Dictionary of receivers at which to calculate arrival times

        """

        available = [receiver
                     for uid, receiver in receivers.items()
                     if receiver.available(self.otime)]

        for receiver in available:
            arrival = Arrival(self, receiver)
            arrival.get_iris_picks(phases, model)
            self + arrival

    def load_waveforms(self):
        """
        Load the waveform data for all Arrivals associated with this Source

        """

        for uid, arrival in self.arrivals.items():
            arrival.load_waveform()

    @property
    def output(self):
        """
        Get a dictionary representation of the Source object for output to
        JSON

        Returns
        -------
        Dictionary containing information about the Source, related Arrivals
        and their Picks in a JSON serialisable format

        """

        arrivals = {}
        for uid, arrival in self.arrivals.items():
            arrivals[uid] = arrival.output

        out = {"otime": str(self.otime),
               "latitude": self.latitude,
               "longitude": self.longitude,
               "depth": self.depth,
               "magnitude": self.magnitude,
               "uid": self.uid,
               "arrivals": arrivals}

        return out

    def __add__(self, arrival):
        """
        Rich operator to associate an Arrival with the Source

        Parameters
        ----------
        arrival : Arrival object
            Contains information about the arrival

        """

        self.arrivals[arrival.receiver.uid] = arrival

    def __sub__(self, arrival):
        """
        Rich operator to remove an associated Arrival from the Source

        Parameters
        ----------
        arrival : Arrival object
            Contains information about the arrival

        """

        del self.arrivals[arrival.receiver.uid]

    def __len__(self):
        """
        Returns the number of Arrivals associated with this Source

        """

        return len(self.arrivals)

    def __iter__(self):
        return iter(self.arrivals.items())
