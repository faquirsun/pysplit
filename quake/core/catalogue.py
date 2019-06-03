# -*- coding: utf-8 -*-
"""
This module requires that "cartopy", "matplotlib" and "obspy" are installed in
your Python environment.

Author: Conor Bacon
"""

import json
import pathlib

import cartopy.crs as ccrs
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel

from quake.core.source import Source


class Catalogue(object):
    """
    Catalogue class

    The Catalogue class provides a means of collecting information on a number
    of seismic Sources. Can be associated with a specific Network and Archive
    and methods are provided to retrieve arrival information and waveform data
    from these, respectively.

    Employs routines from obspy.

    Attributes
    ----------

    Methods
    -------

    """

    def __init__(self, path, archive=None, network=None, new=False):
        """
        Class initialisation

        Parameters
        ----------
        new : bool, optional
            If True, will also generate a metafile (default: False)
        kwargs : dict
            Dictionary containing Catalogue attributes and values

        """

        super().__init__()

        self.path = pathlib.Path(path)

        self.archive = archive
        self.network = network

        if new:
            self._create_directories()
            metapath = self.path / "meta" / self.path.stem
            if self.archive is not None:
                self.archive.save(metapath)
            if self.network is not None:
                self.network.save(metapath)

        self.sources = {}

    def set_params(self, params):

        for key, value in params.items():
            setattr(self, key, value)

        metafile = (self.path / "meta" / self.uid).with_suffix(".params")
        if not metafile.is_file():
            with metafile.open("w") as f:
                json.dump(params, f, sort_keys=True,
                          indent=4, ensure_ascii=False)

    def sources_from_datacentre(self):
        """
        Generates a Catalogue of Sources from the IRIS webclient API.

        """

        client = Client("IRIS")

        st = UTCDateTime(self.sdate)
        et = UTCDateTime(self.edate)

        sources = client.get_events(starttime=st,
                                    endtime=et,
                                    minmagnitude=self.minmag,
                                    latitude=self.network.centre[1],
                                    longitude=self.network.centre[0],
                                    minradius=self.minrad,
                                    maxradius=self.maxrad)

        for source in sources:
            try:
                otime = source.preferred_origin().time
                lat = source.preferred_origin().latitude
                lon = source.preferred_origin().longitude
                dep = source.preferred_origin().depth / 1000.0
                mag = source.preferred_magnitude().mag

                sourceid = otime.isoformat()
                for char_ in [":", "-", "T", ".", " ", "Z"]:
                    sourceid = sourceid.replace(char_, "")

                src_info = {"otime": otime,
                            "latitude": lat,
                            "longitude": lon,
                            "depth": dep,
                            "magnitude": mag,
                            "uid": sourceid}

                self + Source(self, **src_info)
            except TypeError:
                pass

        self.save()

    def arrivals_from_datacentre(self, phases, vmodel="ak135"):
        # def get_arrivals(self, phases=["P", "S"], vmodel="ak135"):
        """
        Generates predicted arrival times using a whole earth velocity model

        Parameters
        ----------
        phases : list of strings, optional
            List containing the seismic phases for which to calculate
            predicted traveltimes. Defaults to "SKS"

        """
        model = TauPyModel(model=vmodel)

        for sourceid, source in self.sources.items():
            source.get_iris_arrivals(phases, model, self.network.receivers)

        self.save()

    def lookup(self, sourceid):
        """
        Queries the sources for a given source.

        Parameters
        ----------
        sourceid : str
            Source name/unique identifier

        Returns
        -------
        Source object

        """

        if type(sourceid) != str:
            print("You must identify the source with a string.")
            return
        if sourceid in self.sources:
            return self.sources[sourceid]
        else:
            msg = "There is no source with id '{}' in this catalogue."
            print(msg.format(sourceid))

    def save(self, path=None):
        """
        Output contents of catalogue to .json file as well as the Network and
        Archive objects

        Parameters
        ----------
        path : pathlib Path object or str
            Location to save file

        """

        catalogue = self.output

        if path is None:
            path = self.path / "meta" / self.uid
        else:
            path = pathlib.Path(path)
        outfile = path.with_suffix(".catalogue")
        with outfile.open("w") as f:
            json.dump(catalogue, f, sort_keys=True, indent=4,
                      ensure_ascii=False)

    def load(self, path=None):
        """
        Load contents of .json file.

        Parameters
        ----------
        path : pathlib Path object or str
            Location of file to load

        """

        if path is None:
            infile = (self.path / "meta" / self.uid).with_suffix(".catalogue")
        else:
            infile = pathlib.Path(path)

        if infile.suffix == ".catalogue" or infile.suffix == ".json":
            if infile.is_file():
                # Read in source .json
                with infile.open("r") as f:
                    catalogue = json.load(f)
                sources = catalogue["sources"]
                # Parse source files to Source objects
                for uid, source in sources.items():
                    src = Source(source, self.network)
                    src.load()
                    self + src
        else:
            print("Please provide a .json file.")

    def plot(self, map_widget, lims=None):
        """
        Creates a geographic representation of the catalogue.

        Parameters
        ----------
        map_widget : Axes object, optional
            The axes object on which to plot the map. One is generated
            if not provided
        lims : list of floats, optional
            Contains the geograhical boundaries of the map. If not provided,
            suitable limits are calculated from the source locations

        """

        nlon, nlat = self.network.centre
        lons = [s.longitude for u, s in self.sources.items()]
        lats = [s.latitude for u, s in self.sources.items()]

        if lims is None:
            londiff = (max(lons) - min(lons)) * 0.1
            latdiff = (max(lats) - min(lats)) * 0.1
            self.lon0 = min(lons) - londiff
            self.lon1 = max(lons) + londiff
            self.lat0 = min(lats) - latdiff
            self.lat1 = max(lats) + latdiff
        else:
            self.lon0 = lims["lon0"]
            self.lon1 = lims["lon1"]
            self.lat0 = lims["lat0"]
            self.lat1 = lims["lat1"]

        if self.type == "teleseismic":
            proj = ccrs.AzimuthalEquidistant(
                central_longitude=nlon,
                central_latitude=nlat)
        elif self.type == "local":
            proj = ccrs.PlateCarree()

        ax = map_widget.canvas.mapPlot(proj)
        if self.type == "local":
            ax.set_xlabel("Longitude, degrees", fontsize=10)
            ax.set_ylabel("Latitude, degrees", fontsize=10)
            ax.set_extent([self.lon0, self.lon1, self.lat0, self.lat1], proj)
        ax.set_aspect("auto")

        self._plot_sources(ax)
        self._plot_receivers(ax)

    @property
    def output(self):
        """
        Get a dictionary representation of the Catalogue object for output to
        JSON

        Returns
        -------
        Dictionary containing information about the Catalogue, related Sources,
        their Arrivals and related Picks in a JSON serialisable format

        """

        sources = {}
        for uid, source in self.sources.items():
            sources[uid] = source.output

        out = {"sources": sources}

        return out

    def __add__(self, source):
        """
        Rich operator to associate a Source with the Catalogue

        Parameters
        ----------
        source : Source object
            Contains information about the source

        """

        self.sources[source.uid] = source

    def __sub__(self, source):
        """
        Rich operator to remove an associated Source from the Catalogue

        Parameters
        ----------
        source : Source object
            Contains information about the source

        """

        del self.sources[source.uid]

    def __len__(self):
        """
        Returns the number of Sources in the Catalogue

        """

        return len(self.sources)

    def _create_directories(self):
        """
        Generate directory structure for a new Catalogue

        """

        dirs = ["meta", "mseed", "tmp", "plots", "picks", "exports"]
        for dir_ in dirs:
            new_dir = self.path / dir_
            new_dir.mkdir(parents=True, exist_ok=True)

    def _plot_sources(self, ax):
        """
        Plot source locations

        Parameters
        ----------
        ax : matplotlib Axes object
            Axes on which to plot source locations

        """

        tolerance = 10
        for uid, source in self.sources.items():
            ax.scatter(source.longitude, source.latitude, s=70,
                       marker="o", color="white", picker=tolerance, zorder=10,
                       alpha=0.9, label="SOURCE: {}".format(uid),
                       transform=ccrs.Geodetic(), edgecolors="k")

    def _plot_receivers(self, ax):
        """
        Plot receiver locations

        Parameters
        ----------
        ax : matplotlib Axes object
            Axes on which to plot receiver locations

        """

        tolerance = 10
        for uid, receiver in self.network.receivers.items():
            ax.scatter(receiver.longitude, receiver.latitude, s=75,
                       marker="v", color="green", picker=tolerance, zorder=15,
                       alpha=0.9, label="RECEIVER: {}".format(uid),
                       transform=ccrs.Geodetic(), edgecolors="k")
