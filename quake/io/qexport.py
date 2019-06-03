# -*- coding: utf-8 -*-
"""
This module requires that "obspy" is installed in your Python environment.

The functions in this module provide a set of export functions that transform
the data contained within a Quake catalogue and generates the relevant input
files for other pieces of software.

Author: Conor Bacon
"""


def export_nonlinloc(catalogue):
    """
    Format of .nonlinloc files:

    Station | ? | ? | ? | PHASE | POLARITY | YEARMONDAY | HOURMIN | SECMSEC

    File names:

    UID.nonlinloc

    TO-DO
    -----
    Establish and implement a pick priority (manual -> qmigrate -> modelled)

    """

    out_fmt = "{station}\t ? \t ? \t ? \t {phase}\t{polarity} {ymd} {hm} {sms}"
    out_fmt += " GAU {error}  0.00e+00  0.00e+00  0.00e+00  1.00e+00\n"

    export_to = catalogue.path / "exports" / "nonlinloc"
    export_to.mkdir(exist_ok=True)

    for i, source in catalogue.sources.items():

        fname = (export_to / "{}".format(source.uid)).with_suffix(".nonlinloc")
        with fname.open("w") as f:
            for i, arrival in source.arrivals.items():
                for pick in ["P_hyp", "S_hyp"]:
                    try:
                        ttime = source.otime + float(arrival.picks[pick].ttime[0][0])
                    except KeyError:
                        continue
                    polarity = arrival.picks[pick].polarity
                    ymd = "{}{}{}".format(str(ttime.year),
                                          str(ttime.month).zfill(2),
                                          str(ttime.day).zfill(2))
                    hm = "{}{}".format(str(ttime.hour).zfill(2),
                                       str(ttime.minute).zfill(2))
                    sms = "{}.{}".format(str(ttime.second).zfill(2),
                                         str(int(ttime.microsecond / 100)).zfill(4))
                    line = out_fmt.format(station=arrival.receiver.uid,
                                          phase=pick[0],
                                          polarity=polarity,
                                          ymd=ymd,
                                          hm=hm,
                                          sms=sms,
                                          error="1.00e-02")

                    f.write(line)
