# -*- coding: utf-8 -*-
"""
There are no requirements for this module.

Author: Conor Bacon
"""


class Pick(object):
    """
    Pick class

    The Pick class encapsulates the information that uniquely describes a
    seismic phase arrival at a receiver, including the pick identifier code,
    traveltime, pick error and pick polarity.

    Attributes
    ----------
    uid : str
        The pick code (e.g. seismic phase and indication of origin of pick)
    ttime : float
        Traveltime relative to origin time
    error : float
        Error on traveltime pick
    polarity : str
        First-motion polarisation of the phase arrival

    """

    def __init__(self, pickid, ttime, error=None, polarity="?"):
        """
        Class initialisation

        Parameters
        ----------
        pickid : str

        time : float

        method : str


        """

        self.uid = pickid
        if isinstance(ttime, list):
            self.ttime = [] + ttime
        else:
            self.ttime = [ttime]
        self.error = error
        self.polarity = polarity

    @property
    def output(self):
        """
        Get a dictionary representation of the Pick object for output to JSON.

        """

        return vars(self)
