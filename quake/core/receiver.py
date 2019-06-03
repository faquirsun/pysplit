# -*- coding: utf-8 -*-
"""
This module requires that "obspy" is installed in your Python environment.

Author: Conor Bacon
"""

from obspy import UTCDateTime


class Receiver(object):
    """
    Receiver class

    The Receiver class encapsulates the information that uniquely describes a
    receiver, including the network/station identifier codes, geographical
    location and deployment period.

    Attributes
    ----------
    network : str, optional
        The network code (default is an empty string)
    uid : str
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

    Methods
    -------
    available(time)
        Queries the Receiver availability for the specified time period

    """

    def __init__(self, rec_info, network=None):
        """
        Class initialisation

        Parameters
        ----------
        rec_info : pandas DataFrame
            Contains the information to initialise a receiver object

        """

        super().__init__()

        for key, value in rec_info.items():
            setattr(self, key, value)

        if network is not None:
            self.network = network
        elif not hasattr(self, "network"):
            self.network = ""

        if hasattr(self, "deployment"):
            self.deployment = UTCDateTime(self.deployment)
        else:
            self.deployment = ""
        if hasattr(self, "retrieval"):
            self.retrieval = UTCDateTime(self.retrieval)
        else:
            self.retrieval = ""

        if not hasattr(self, "arrivals"):
            self.arrivals = []

    def available(self, time):
        """
        Queries the Receiver availability for the specified time period.

        Parameters
        ----------
        time : UTCDateTime object
            DateTime against which to check availability

        Returns
        -------
        True if available otherwise False

        """

        if isinstance(self.deployment, str) or isinstance(self.retrieval, str):
            return True

        if time.between(self.deployment, self.retrieval):
            return True
        else:
            return False

    @property
    def output(self):
        """
        Get a dictionary representation of the Receiver object for output to
        JSON

        Returns
        -------
        Dictionary containing information about the Receiver and related
        Arrival uids in a JSON serialisable format

        """

        output = {"network": self.network,
                  "uid": self.uid,
                  "longitude": self.longitude,
                  "latitude": self.latitude,
                  "elevation": self.elevation,
                  "deployment": str(self.deployment),
                  "retrieval": str(self.retrieval),
                  "arrivals": self.arrivals}

        return output

    def __add__(self, source):
        """
        Rich operator to associate an Arrival with the Receiver.

        Parameters
        ----------
        source : Source object
            Contains information about the source

        """

        self.arrivals.append(source.uid)

    def __sub__(self, source):
        """
        Rich operator to remove an associated Arrival from the Receiver.

        Parameters
        ----------
        source : Source object
            Contains information about the source

        """
        if source.uid in self.arrivals:
            self.arrivals.remove(source.uid)
