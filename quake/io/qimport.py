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
import pandas as pd

from quake.core.arrival import Arrival
from quake.core.pick import Pick
from quake.core.source import Source


def parse(input_, catalogue, network, source=None):
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
                qm = input_.glob("*.event")
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

        picks = pd.read_csv(source.with_suffix(".stn"))
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
