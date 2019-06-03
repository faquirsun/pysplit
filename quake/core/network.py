# -*- coding: utf-8 -*-
"""
This module requires that "pandas" is installed in your Python environment.

Author: Conor Bacon
"""

import json
import pathlib

import pandas as pd

from quake.core.receiver import Receiver


class Network(object):
    """
    Network class

    The Network class encapsulates the information that uniquely describes a
    network, including the network identifier code, a geographical centre and
    a composition of Receiver objects associated with the Network.

    Attributes
    ----------
    uid : str, optional
        The network code (default is an empty string)
    receivers : dict
        Dictionary of receivers in the network, with the station UID as keys
    centre : [float, float]
        Geographical location corresponding to the centre of the network

    Methods
    -------
    lookup(receiverid)
        Queries the network with the identifier and returns a Receiver object

    """

    def __init__(self, receiver_file=None, network_uid=""):
        """
        Class initialisation

        Parameters
        ----------
        receiver_file : pathlib Path object or str
            Path pointing to json or csv file containing receiver information

        """

        super().__init__()

        self.uid = network_uid
        self.receivers = {}

        if receiver_file is not None:
            self.load(receiver_file)

    def lookup(self, receiverid):
        """
        Queries the network for a given receiver

        Parameters
        ----------
        receiverid : str
            Receiver name/unique identifier

        Returns
        -------
        Receiver object

        """

        if not type(receiverid) == str:
            print("You must identify the receiver with a string.")
            return
        if receiverid in self.receivers:
            return self.receivers[receiverid]
        else:
            msg = "There is no receiver with id '{}' in this network."
            print(msg.format(receiverid))

    def save(self, path):
        """
        Output contents of network to .json file

        Parameters
        ----------
        path : pathlib Path object or str
            Location to save file

        """

        out = self.output

        path = pathlib.Path(path)
        outfile = path.with_suffix(".network")
        with outfile.open("w") as f:
            json.dump(out, f, sort_keys=True, indent=4, ensure_ascii=False)

    def load(self, path):
        """
        Load contents of .json/.csv file

        Parameters
        ----------
        path : pathlib Path object or str
            Location of file to load

        """

        infile = pathlib.Path(path)
        if infile.suffix == ".network" or infile.suffix == ".json":
            with infile.open("r") as f:
                network = json.load(f)
            self.uid = network["uid"]
            for uid, receiver in network["receivers"].items():
                receiver = Receiver(receiver)
                self + receiver
            self.centre = network["centre"]
        else:
            receivers = pd.read_csv(infile)
            for receiver in receivers.to_dict("records"):
                receiver = Receiver(receiver)
                self + receiver
            self.centre = self._get_centre()

    @property
    def output(self):
        """
        Get a dictionary representation of the Network object for output to
        JSON

        Returns
        -------
        Dictionary containing information about the Network and related
        Receivers in a JSON serialisable format

        """

        receivers = {}
        for uid, receiver in self.receivers.items():
            receivers[uid] = receiver.output

        out = {"uid": self.uid,
               "centre": self.centre,
               "receivers": receivers}

        return out

    def __add__(self, receiver):
        """
        Rich operator to associate a Receiver with the Network

        Parameters
        ----------
        receiver : Receiver object
            Contains information about the receiver

        """

        self.receivers[receiver.uid] = receiver

    def __sub__(self, receiver):
        """
        Rich operator to remove an associated Receiver from the Network

        Parameters
        ----------
        receiver : Receiver object
            Contains information about the receiver

        """

        del self.receivers[receiver.uid]

    def __len__(self):
        """
        Returns the number of Receivers in the Network

        """

        return len(self.receivers)

    def _get_centre(self):
        """Get the lon/lat centre of the network"""
        lons = [s.longitude for u, s in self.receivers.items()]
        lats = [s.latitude for u, s in self.receivers.items()]

        clon = (max(lons) + min(lons)) / 2
        clat = (max(lats) + min(lats)) / 2

        return [clon, clat]
