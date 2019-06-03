# -*- coding: utf-8 -*-
"""
This module has no requirements.

Author: Conor Bacon
"""

import json
import pathlib


class Archive(object):
    """
    Archive object

    The Archive class encapsulates the information required to specify the
    format and location of locally stored mSEED data.

    Attributes
    ----------
    path : pathlib Path object
        Location of raw mSEED files
    schema : str
        Identifying code for archive format schema
    structure : str
        File naming format within archive

    Methods
    -------
    archive_structure(schema, custom)
        Set the file naming format of the data archive
    save(path)

    """

    def __init__(self, path):
        """
        MSEED object initialisation

        Parameters
        ----------
        path : str
            Location of raw mSEED files

        """

        self.path = pathlib.Path(path)

    def archive_structure(self, schema="SeisComp3", custom=None):
        """
        Define the format of the data archive

        Parameters
        ----------
        schema : str, optional
            Sets path type for different archive formats
        custom : str
            Custom formattable string that represents how data is stored in the
            archive

        """

        self.schema = schema

        if schema == "SeisComp3":
            self.structure = "{year}/*/{receiver}/*/*.{receiver}..*.D.{year}.{jday}"
        elif schema == "Cambridge Volcano Seismology":
            self.structure = "{year}/{jday}/*_{receiver}_*"
        elif schema == "QuakeMigrate":
            self.structure = "*_{sourceid}"
        elif schema == "YEAR/JD/STATION":
            self.structure = "{year}/{jday}/{receiver}*"
        elif schema == "STATION.YEAR.JULIANDAY":
            self.structure = "*{receiver}.*.{year}.{jday}"
        elif schema == "/STATION/STATION.YearMonthDay":
            self.structure = "{receiver}/{receiver}.{year}{month:02d}{day:02d}"
        elif schema == "YEAR_JD/STATION":
            self.structure = "{year}_{jday}/{receiver}_*"
        elif schema == "custom":
            self.structure = custom

    def save(self, path):
        """
        Output contents of Archive to .archive JSON file

        Parameters
        ----------
        path : pathlib Path object or str
            Location to save file

        """

        out = self.output

        path = pathlib.Path(path)
        outfile = path.with_suffix(".archive")
        with outfile.open("w") as f:
            json.dump(out, f, sort_keys=True, indent=4, ensure_ascii=False)

    def load(self, path):
        """
        Load contents of .archive JSON file

        Parameters
        ----------
        path : pathlib Path object or str
            Location of file to load

        """

        infile = pathlib.Path(path)
        if infile.suffix == ".archive":
            with infile.open("r") as f:
                archive = json.load(f)
            self.path = pathlib.Path(archive["path"])
            self.schema = archive["schema"]
            self.structure = archive["structure"]

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

        out = {"path": str(self.path),
               "schema": self.schema,
               "structure": self.structure}

        return out
