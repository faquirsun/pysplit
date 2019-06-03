# -*- coding: utf-8 -*-
"""
This module requires that "obspy" and "numpy" are installed in your Python
environment.

Author: Conor Bacon
"""

from itertools import chain

import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import locations2degrees
from obspy import read, Stream

from quake.core.pick import Pick


class Arrival(object):
    """
    Arrival class

    The Arrival class encapsulates the information that uniquely describes the
    arrival of a seismic signal at a specific receiver.

    Records the component data as well as the pick times of phases (both
    automated and manually refined).

    It forms an aggregation of the information from a particular Source object
    and a particular Receiver object.

    Attributes
    ----------
    file_format : str
        formattable string used for naming output files
    source : Source-type object
        a Source object that contains information specific to the source
    receiver : Receiver-type object
        a Receiver object that contains information specific to the receiver
    components : str
        a string that tracks the current component system ("ZNE" or "ZRT")
    pick_path : str
        a formatted string that points to the location of the relevant pick
        files

    Methods
    -------
    filter(code=None, filt)
        Filters the seismic data
    rotate(method)
        Rotates the seismic data
    saveData()
        Saves manual pick data
    addData(data, value, pick_type=None)
        Adds pick information to dictionary

    Properties
    ----------
    network : string, optional
        Network code (default is an empty string)
    station : string, optional
        Station code
    location : string, optional
        Location code (default: "")
    starttime : UTCDateTime, optional
        Date and time of the first data sample given in UTC
    endtime : UTCDateTime, optional
        Date and time of last data sample given in UTC
    sampling_rate : float, optional
        Sampling rate (units: Hz; default value: 1.0)
    delta : float, optional
        Sample distance (units: s; default value: 1.0)
    npts : int, optional
        Number of sample points (default value: 0)
    calib : float, optional
        Calibration factor (default value: 1.0)

    """

    uid_stub = "{}_{}"

    def __init__(self, source, receiver, load_stream=False):
        """
        Object initialisation

        Parameters
        ----------
        catalogue_path : str
            Path to catalogue containing this source-receiver pair
        source : Source object
            Contains information specific to the source
        receiver : Receiver object
            Contains information specific to the receiver
        """

        self.source = source
        self.receiver = receiver

        self.uid = self.uid_stub.format(source.uid, receiver.uid)

        # self.path = source.catalogue.path

        self.components = "ZNE"

        if load_stream:
            self.load_stream()

        self.picks = {}

        self.distance, self.az, self.baz = gps2dist_azimuth(source.latitude,
                                                            source.longitude,
                                                            receiver.latitude,
                                                            receiver.longitude)

    def load_picks(self, picks):
        """

        """

        for uid, pick in picks.items():
            pk = Pick(uid, pick["ttime"], pick["error"], pick["polarity"])
            self + pk

    def get_iris_picks(self, phases, model):
        dist = locations2degrees(
            self.source.latitude, self.source.longitude,
            self.receiver.latitude, self.receiver.longitude)

        phases = model.get_travel_times(
            source_depth_in_km=self.depth,
            distance_in_degree=dist,
            phase_list=phases,
            receiver_depth_in_km=self.receiver.elevation / 1000)

        if not phases:
            return
        else:
            for phase in phases:
                pick_id = "{}_model".format(phase.purist_name)
                time = phase.time

                pick = Pick(pick_id, time)

                if pick_id in self.picks:
                    self.picks[pick_id].ttime.append(time)
                else:
                    self + pick

    # def load_stream(self):
    #     """
    #     Load streams

    #     """
    #     path = self.path / "data" / self.receiver.uid / self.uid
    #     path = path.with_suffix(".*")

    #     self.stream = read(str(path))
    #     self._update_components(self.stream)

    @property
    def output(self):
        """
        Get a dictionary representation of the Source object for output to
        JSON

        """

        picks = {}
        for uid, pick in self.picks.items():
            picks[uid] = pick.output

        out = {"sourceid": self.source.uid,
               "receiverid": self.receiver.uid,
               "picks": picks}

        return out

    @property
    def max_ttime(self):
        ttimes = []
        for uid, pick in self.picks.items():
            for ttime in pick.ttime:
                ttimes += [float(ttime)]

        if not ttimes:
            window = 120.
        else:
            window = max(ttimes)

        return window

    def __add__(self, pick):
        """
        Rich operator to associate a Pick with the Arrival.

        Parameters
        ----------
        pick : Pick object
            Contains information about the pick

        """

        self.picks[pick.uid] = pick

    def __sub__(self, pick):
        """
        Rich operator to remove an associated Pick from the Arrival

        Parameters
        ----------
        pick : Pick object
            Contains information about the pick

        """

        del self.picks[pick.uid]

    def __iter__(self):
        return iter(self.picks.items())

    def delete(self):
        path = self.path / "data" / self.receiver.uid
        files = path.glob("{}.*".format(self.uid))
        for file in files:
            file.unlink()

        self.source.delete_arrival(self.receiver.uid)

    def filter(self, method="obspy", filt=None):
        """
        Filters the seismic data and updates the component attributes

        Parameters
        ----------
        method : str, optional
            The method to be used for filtering. Options are "obspy", "sac" or
            None (default to None, which removes the filter)
        filt : dict
            Contains information about the filter

        """

        tmp_stream = self.stream.copy()

        if filt is None:
            (tmp_stream.detrend("linear")).detrend("demean")

        else:
            if method == "sac":
                pass
                # if filt["zerophase"]:
                #     ph = 2
                # else:
                #     ph = 1
                # sac_cmd = "sac << EOF >/dev/null\n"
                # sac_cmd += "r f.sac\n"
                # sac_cmd += "rmean\n"
                # sac_cmd += "taper\n"
                # sac_cmd += "bp c {:f} {:f} n {:d} p {:d}\n"
                # sac_cmd += "write over\n"
                # sac_cmd += "q\n"
                # sac_cmd += "EOF\n"

                # tmp_dir = pathlib.Path(self.path / "tmp")
                # st = Stream()
                # for tr in tmp_stream:
                #     tr.write("{}/f.sac".format(tmp_dir), "SAC")
                #     text = sac_cmd.format(filt["freqmin"], filt["freqmax"],
                #                           filt["corners"], ph)
                #     os.chdir(tmp_dir)
                #     open("s", "w+").write(text)
                #     os.system("bash s")
                #     tr = read("{}/f.sac".format(tmp_dir))[0]
                #     st += tr
            elif method == "obspy":
                tmp_stream.taper(max_percentage=0.05)
                tmp_stream.filter(**filt)

        self._update_components(tmp_stream)

    def rotate(self, method=None, phi=None):
        """
        Rotates the seismic data and updates the component attributes

        Parameters
        ----------
        method : str
            Specifies the rotation operation to be used
        phi : float
            Angle through which to rotate components (units: degrees)

        """

        if method is None:
            rphi = np.deg2rad(phi)
            # Finish arbitrary rotation

        elif method == "NE->RT":
            self.components = "ZRT"
        elif method == "RT->NE":
            self.components = "ZNE"

        self.stream.rotate(method, back_azimuth=self.baz)

        self._update_components(self.stream)

    def plot(self, axes=None, lims=None):
        """
        Plot component traces on provided axes with specified limits

        Parameters
        ----------
        axes : list of Axes objects, optional
            An array containing the axes
        lims : list of float lists

        """

        if lims is None:
            xlims = (self.comp1.times(type="relative")[0],
                     self.comp1.times(type="relative")[-1])
            ylims = (-1, 1)

        else:
            xlims = lims[0]
            ylims = lims[1]

        for ax in axes:
            ax.set_xlim(xlims)
            ax.set_ylim(ylims)

        axes[0].plot(self.comp1.times(type="relative"),
                     self.comp1.data / self.norm_factor,
                     linewidth=0.5, c="white")
        axes[1].plot(self.comp2.times(type="relative"),
                     self.comp2.data / self.norm_factor,
                     linewidth=0.5, c="white")
        axes[2].plot(self.comp3.times(type="relative"),
                     self.comp3.data / self.norm_factor,
                     linewidth=0.5, c="white")

        axes[0].set_ylabel("Z")
        if self.components == "ZNE":
            axes[1].set_ylabel("N")
            axes[2].set_ylabel("E")

        elif self.components == "ZRT":
            axes[1].set_ylabel("R")
            axes[2].set_ylabel("T")

    def spectrogram(self, ax=None, log=True):
        """
        Plot spectrogram of the stream data

        Parameters
        ----------
        ax : Axes object, optional
            Axes object on which the spectrogram is plotted
        log : Boolean
            Boolean variable that specifies whether to use a log scale
            on the y-axis of the spectrogram (default True)

        """

        self.stream.spectrogram(log=log, axes=ax)
        ax.set_ylabel("log(Power)")

    def hodogram(self, ax, phi=None, dt=None):
        if not phi and not dt:
            x = self.comp3.data
            y = self.comp2.data
        else:
            x, y = self.correctSplitting(phi, dt)

        # Normalise the particle motion plots
        norm = max(max(abs(x)), max(abs(y)))

        x = [z / norm for z in x]
        y = [z / norm for z in y]

        # Plot to axes
        ax.plot(x, y)

    @property
    def components(self):
        """Get and set the component coordinate system"""
        return self._components

    @components.setter
    def components(self, value):
        self._components = value

    @property
    def distance(self):
        """Get and set the epicentral distance between source and receiver"""
        return self._distance

    @distance.setter
    def distance(self, value):
        self._distance = value

    @property
    def az(self):
        """Get and set the azimuth between source and receiver"""
        return self._az

    @az.setter
    def az(self, value):
        self._az = value

    @property
    def baz(self):
        """Get and set the back-azimuth between source and receiver"""
        return self._baz

    @baz.setter
    def baz(self, value):
        self._baz = value

    @property
    def norm_factor(self):
        """Get the normalising factor for the traces"""
        z_max = max(abs(self.comp1.data))
        n_max = max(abs(self.comp2.data))
        e_max = max(abs(self.comp3.data))

        return max(z_max, n_max, e_max)

    @property
    def network(self):
        """Get the network code"""
        return self.receiver.network

    @property
    def station(self):
        """Get the station code"""
        return self.receiver.station

    @property
    def location(self):
        """Get the geographical coordinates of the source and receiver"""
        return self.source.longitude, self.source.latitude, \
            self.source.depth, self.receiver.lon, self.receiver.lat, \
            self.receiver.elv

    @property
    def starttime(self):
        """Get the date and time of the first data sample"""
        return self.comp1.stats.starttime

    @property
    def endtime(self):
        """Get the date and time of the last data sample"""
        return self.comp1.stats.endtime

    @property
    def sampling_rate(self):
        """Get the sampling rate"""
        return self.comp1.stats.sampling_rate

    @property
    def delta(self):
        """Get the sample distance"""
        return self.comp1.stats.delta

    @property
    def npts(self):
        """Get the number of data"""
        return self.comp1.stats.npts

    @property
    def calib(self):
        """Get the calibration factor"""
        return self.comp1.stats.calib

    def load_waveform(self, archive, outpath):
        """
        Retrieves waveform data for a single source-receiver pair from an
        archive.

        Parameters
        ----------

        """

        # Check first to see if already retrieved
        tmp = outpath / "mseed" / self.receiver.uid / self.source.uid
        print(tmp)
        if tmp.with_suffix(".e").is_file():
            print("Hello")
            st = read(str(tmp.with_suffix(".e")))
            st += read(str(tmp.with_suffix(".n")))
            st += read(str(tmp.with_suffix(".z")))

            self.stream = st

            self._update_components(self.stream)
            return

        # Get max traveltime by searching through available picks. If none,
        # simply take base window of 300 seconds?
        window = self.max_ttime

        wbeg = self.source.otime - 60.
        wend = self.source.otime + window

        events = archive.structure.format(year=wbeg.year,
                                          month=wbeg.month,
                                          jday=str(wbeg.julday).zfill(3),
                                          receiver=self.receiver.uid,
                                          comp="*")
        files = archive.path.glob(events)
        if wbeg.julday != wend.julday:
            events = archive.structure.format(year=wend.year,
                                              month=wend.month,
                                              jday=str(wend.julday).zfill(3),
                                              receiver=self.receiver.uid,
                                              comp="*")

            files = chain(files, archive.path.glob(events))

        st = Stream()
        try:
            # "Peep" the iterator for available files
            first = next(files)
            files = chain([first], files)
            for file in files:
                file = str(file)
                try:
                    st += read(file)
                except TypeError:
                    msg = "File not in mSEED format."
                    print(msg)
                    continue

            # Do some cleanup
            st.merge(method=-1)
            st.trim(starttime=wbeg, endtime=wend)

            if len(st) != 3:
                print("Archive incomplete for this arrival.")
                # Archive incomplete for this time period - do not write
            else:
                # Write to files
                for tr in st:
                    name = "{stat}/{sourceid}.{comp}"
                    out = outpath / "mseed" / name.format(
                        stat=self.receiver.uid,
                        sourceid=self.source.uid,
                        comp=tr.stats.channel[-1].lower())

                    tr.write(str(out), format="MSEED")

                self.stream = st

        except StopIteration:
            print("No files available in the archive for this arrival.")
            # Do some more handling here?

    def _update_components(self, stream):
        """Update the component attributes

        Parameters
        ----------
        stream : Stream object
            List like object of multiple ObsPy Trace objects

        """

        if self.components == "ZNE":
            self.comp1 = stream.select(channel="*Z")[0]
            self.comp2 = stream.select(channel="*N")[0]
            self.comp3 = stream.select(channel="*E")[0]
        elif self.components == "ZRT":
            self.comp1 = stream.select(channel="*Z")[0]
            self.comp2 = stream.select(channel="*R")[0]
            self.comp3 = stream.select(channel="*T")[0]
