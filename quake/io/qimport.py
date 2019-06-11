# -*- coding: utf-8 -*-
"""
This module requires that "obspy" is installed in your Python environment.

The functions in this module provide a set of import functions that parse the
data from various file formats into a Quake catalogue.

Author: Conor Bacon
"""

from itertools import chain
import pathlib

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import locations2degrees
from obspy.taup import TauPyModel
import pandas as pd

from quake.core.arrival import Arrival
from quake.core.pick import Pick
from quake.core.source import Source


def parse_local(input_, catalogue, network, source=None):
    """
    If source is None, tries to guess the data format that the user is trying
    to parse and redirects to the correct parser

    """

    input_ = pathlib.Path(input_)
    if source is not None:
        if source == ".hyp":
            _import_hyp(input_, catalogue, network)
        elif source == ".nonlinloc":
            _import_nll(input_, catalogue, network)
        elif source == "QuakeMigrate":
            _import_qmigrate(input_, catalogue, network)
        elif source == "IRIS":
            _import_iris(input_, catalogue, network)
    else:
        if input_.is_file():
            if input_.suffix == ".hyp":
                _import_hyp(input_, catalogue, network)
            else:
                msg = "{} file format is not currently supported."
                print(msg.format(input_.suffix))
        elif input_.is_dir():
            # Parser for QuakeMigrate
            try:
                qm = input_.glob("events/*.event")
                first = next(qm)
                qm = chain([first], qm)
                _import_qmigrate(qm, catalogue, network)
            except StopIteration:
                pass

            # Parser for folder containing .nonlinloc files
            try:
                nll = input_.glob("*.nonlinloc")
                first = next(nll)
                nll = chain([first], nll)
                _import_nll(nll, catalogue, network)
                return
            except StopIteration:
                pass


def parse_teleseismic(catalogue, network):
    """

    """

    _teleseismic_sources(catalogue, network)
    _teleseismic_arrivals(catalogue, network, phases=["P", "S", "SKS"])


def _teleseismic_sources(catalogue, network):
    """
    Generates a Catalogue of Sources from the IRIS webclient API.

    """

    client = Client("IRIS")

    st = UTCDateTime(catalogue.sdate)
    et = UTCDateTime(catalogue.edate)

    sources = client.get_events(starttime=st,
                                endtime=et,
                                minmagnitude=catalogue.minmag,
                                latitude=catalogue.network.centre[1],
                                longitude=catalogue.network.centre[0],
                                minradius=catalogue.minrad,
                                maxradius=catalogue.maxrad)

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

            catalogue + Source(src_info, network)
        except TypeError:
            pass

    catalogue.save()


def _teleseismic_arrivals(catalogue, network, phases, vmodel="ak135"):
    """
    Generates predicted arrival times using a whole earth velocity model

    Parameters
    ----------
    phases : list of strings, optional
        List containing the seismic phases for which to calculate
        predicted traveltimes. Defaults to "SKS"

    """
    model = TauPyModel(model=vmodel)

    for sourceid, source in catalogue.sources.items():
        _get_iris_arrivals(source, phases, model, network.receivers)

    catalogue.save()


def _get_iris_arrivals(source, phases, model, receivers):
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
                 if receiver.available(source.otime)]

    for receiver in available:
        arrival = Arrival(source, receiver)
        _get_iris_picks(arrival, phases, model)
        source + arrival


def _get_iris_picks(arrival, phases, model):
    dist = locations2degrees(
        arrival.source.latitude, arrival.source.longitude,
        arrival.receiver.latitude, arrival.receiver.longitude)

    phases = model.get_travel_times(
        source_depth_in_km=arrival.source.depth,
        distance_in_degree=dist,
        phase_list=phases,
        receiver_depth_in_km=arrival.receiver.elevation / 1000)

    if not phases:
        return
    else:
        for phase in phases:
            print(phase)
            pick_id = "{}_model".format(phase.purist_name)
            time = phase.time

            pick = Pick(pick_id, time)

            if pick_id in arrival.picks:
                arrival.picks[pick_id].ttime.append(time)
            else:
                arrival + pick


def _import_hyp(hyp, catalogue, network):

    sources = {}
    with hyp.open(mode="r") as f:
        for line in f:
            if "NLLOC " in line:
                picks = []
            elif "GEOGRAPHIC" in line:
                sourceline = line
            elif " > " in line and "PHASE ID" not in line:
                picks.append(line)
            if "END_NLLOC" in line:
                sources[sourceline] = picks

        for src, picks in sources.items():
            src = src.rstrip().split()

            yr, mnth, dy = int(src[2]), int(src[3]), int(src[4])
            hr, mnte, sc = int(src[5]), int(src[6]), float(src[7])
            sec = int(sc)
            msec = int(1e6 * (sc - sec))
            otime = UTCDateTime(yr, mnth, dy, hr, mnte, sec, msec)

            lat = float(src[9])
            lon = float(src[11])
            dep = float(src[13])
            mag = "?"
            sourceid = otime.isoformat()
            for char_ in [":", "-", "T", ".", " "]:
                sourceid = sourceid.replace(char_, "")

            src_info = {"otime": otime,
                        "latitude": lat,
                        "longitude": lon,
                        "depth": dep,
                        "magnitude": mag,
                        "uid": sourceid}

            source = Source(src_info, network)
            # Check if a Source object with same UID is in catalogue
            if source.uid in catalogue.sources:
                source = catalogue.sources[source.uid]
            else:
                catalogue + source

            for pick in picks:
                pick = pick.rstrip().split()
                receiver = network.lookup(pick[0])
                if receiver is None:
                    print("{} not in network.".format(pick[0]))
                    continue

                # Add pick
                pick_id = "{}_hyp".format(pick[4])
                ttime = pick[15]
                pick = Pick(pick_id, ttime)

                arrival = Arrival(source, receiver)
                # Check if an Arrival object with same UID is in catalogue
                if receiver.uid in source.arrivals:
                    arrival = source.arrivals[receiver.uid]
                else:
                    source + arrival

                arrival + pick

    catalogue.save()


def _import_nll(nll_files, catalogue, network):
    pass


def _import_qmigrate(sources, catalogue, network):
    """

    """

    for source in sources:
        otime, _, lon, lat, dep, *_ = pd.read_csv(source).iloc[0]
        otime = UTCDateTime(otime)
        mag = ""
        sourceid = otime.isoformat()
        for char_ in [":", "-", "T", ".", " ", "Z"]:
            sourceid = sourceid.replace(char_, "")

        src_info = {"otime": otime,
                    "latitude": lat,
                    "longitude": lon,
                    "depth": dep,
                    "magnitude": mag,
                    "uid": sourceid}

        src = Source(src_info, network)
        # Check if a Source object with same UID is in catalogue
        if src.uid in catalogue.sources:
            src = catalogue.sources[src.uid]
        else:
            catalogue + src

        pick_file = source.parents[1] / "picks" / "{}".format(source.stem)
        picks = pd.read_csv(pick_file.with_suffix(".picks"))
        for i, pick in picks.iterrows():
            receiver = network.lookup(pick["Name"])
            if receiver is None:
                print("{} not in network.".format(pick["Name"]))
                continue

            arrival = Arrival(src, receiver)
            # Check if an Arrival object with same UID is in catalogue
            if receiver.uid in src.arrivals:
                arrival = src.arrivals[receiver.uid]
            else:
                src + arrival

            # Add picks
            model_id = "{}_model".format(pick["Phase"])
            model_ttime = UTCDateTime(pick["ModelledTime"]) - otime
            model_pick = Pick(model_id, model_ttime)
            arrival + model_pick

            if pick["PickTime"] != "-1":
                qm_id = "{}_qmigrate".format(pick["Phase"])
                qm_ttime = UTCDateTime(pick["PickTime"]) - otime
                pick_error = pick["PickError"]
                qm_pick = Pick(qm_id, qm_ttime, error=pick_error)
                arrival + qm_pick

    catalogue.save()


def _import_iris(iris_file, catalogue, network):
    pass
