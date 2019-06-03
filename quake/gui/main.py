# -*- coding: utf-8 -*-
"""
This module requires that "numpy", "pandas", "matplotlib" and "PyQT5" are
installed in your Python environment

The back-end is powered by a series of classes in quake.core

Author: Conor Bacon
"""

import json
import pathlib
import sys

from matplotlib.patches import Rectangle
import numpy as np
from PyQt5.QtCore import QDateTime, QFile, Qt, QTextStream
import PyQt5.QtWidgets as qt
import PyQt5.QtGui as QtGui
from PyQt5 import uic

from quake.core.archive import Archive
from quake.core.catalogue import Catalogue
from quake.core.network import Network
from quake.gui.picker import Picker
import quake.io.qimport as qimport
import quake.io.qexport as qexport


class Quake(qt.QMainWindow):
    """
    Quake class

    Base class that connects the PyQt5 Quake GUI to the back-end system of
    classes

    Attributes
    ----------

    Methods
    -------

    """

    def __init__(self):
        """
        Class initialisation

        """

        super().__init__()

        self.default_save = False
        self.filt = None

        self.initUI()

    def initUI(self):
        """
        Load UI from designer file and initialise interface

        """

        uic.loadUi("gui/ui/main.ui", self)

        self.connect()

        self.uiMapMpl.hide()

        self.setWindowTitle("Quake - catalogue")
        self.setWindowIcon(QtGui.QIcon("gui/ui/resources/icon.png"))
        self.show()

    def connect(self):
        """
        Connect all signals from Qt to functions

        """

        self.uiNewCatalogueAction.triggered.connect(self.newCatalogue)
        self.uiLoadCatalogueAction.triggered.connect(self.loadCatalogue)
        self.uiDefaultFilterAction.triggered.connect(self.defaultFilter)

        self.uiExportNLL.triggered.connect(lambda: self.export("NLLoc"))
        self.uiExportHypoDD.triggered.connect(lambda: self.export("HypoDD"))
        self.uiExportFMTomo.triggered.connect(lambda: self.export("FMTomo"))

        self.uiMapMpl.canvas.mpl_connect("pick_event",
                                         self._onPick)
        self.uiMapMpl.canvas.mpl_connect("button_press_event",
                                         self._onMapClick)
        self.uiMapMpl.canvas.mpl_connect("motion_notify_event",
                                         self._onMapMove)
        self.uiMapMpl.canvas.mpl_connect("button_release_event",
                                         self._onMapRelease)

        self.uiSourcesListView.clicked.connect(self.sourceSelect)
        self.uiSourcesListView.doubleClicked.connect(self.sourcePick)
        self.uiSourceLoadButton.clicked.connect(self.loadSources)
        self.uiPlotMapButton.clicked.connect(
            lambda: self.plotCatalogueMap(replot=True))
        self.uiResetMapButton.clicked.connect(self.plotCatalogueMap)

    def defaultFilter(self):
        """
        Open a new dialogue window to set a default filter

        If dialogue is not successfully executed, no default filter is set.

        """

        dlg = DefaultFilterDialogue(self)

        if not dlg.exec_():
            return

        self.filt = dlg.filt

    def newCatalogue(self):
        """
        Opens a new dialogue window to generate a new catalogue

        If dialogue is not successfully executed, no catalogue is generated.

        """

        dlg = NewCatalogueDialogue()

        if not dlg.exec_():
            return

        path = dlg.path
        params = dlg.params
        archive = dlg.archive
        network = Network(receiver_file=dlg.network)

        self.createCatalogue(path=path, params=params, archive=archive,
                             network=network, new=True)

    def loadCatalogue(self):
        """
        Opens a QFileDialog window and loads an existing catalogue.

        A test is established to ensure the selected directory is compatible
        with Quake.

        """

        msg = "Choose catalogue directory"
        path = qt.QFileDialog.getExistingDirectory(self, msg)
        path = pathlib.Path(path)

        uid = path.stem

        file = path / "meta" / uid
        try:
            self.uiStatusBar.showMessage("Parsing catalogue parameter file...")
            with file.with_suffix(".params").open("r") as f:
                params = json.load(f)
            network = Network()
            network.load(file.with_suffix(".network"))
            self.uiStatusBar.showMessage("Parsing archive file...")
            with file.with_suffix(".archive").open("r") as f:
                archive = json.load(f)

            self.createCatalogue(path=path, params=params, network=network,
                                 archive=archive)
        except FileNotFoundError:
            msg = "Some of the necessary input .json files are missing."
            qt.QMessageBox.about(self, "Error!", msg)
            return

    def createCatalogue(self, path, archive, network, params, new=False):
        """
        Creates a Catalogue object from the catalogue parameter dictionary

        Parameters
        ----------
        path : str
            Path to where the catalogue directory is set up
        params : dict
            Contains the parameters that define the catalogue
        new : bool, optional
            If True, a new catalogue will be generated, including generation
            of all relevant directories and files (default: False)

        """

        # Interface updates
        self.uiCatalogueNameDisplay.setText(params["uid"])
        self.uiSourcesListView.setEnabled(True)
        if params["type"] == "teleseismic":
            self.uiPlotOptionsStacked.setCurrentIndex(1)
        elif params["type"] == "local":
            self.uiPlotOptionsStacked.setCurrentIndex(0)

        self.archive = Archive(archive["path"])
        self.archive.archive_structure(schema=archive["schema"],
                                       custom=archive["structure"])

        self.network = network

        msg = "Creating catalogue..."
        self.uiStatusBar.showMessage(msg)
        self.catalogue = Catalogue(path, self.archive, network, new)

        self.catalogue.set_params(params)

        self.loadSources()

        if params["type"] == "teleseismic":
            self._populateCatalogueInformation(self.catalogue)

        self.plotCatalogueMap()

    def loadSources(self):
        """
        Loads an existing file or generates one

        """

        self.catalogue.load()

        if len(self.catalogue) == 0:
            # Attempt to generate sources
            msg = "Generating sources from file..."
            self.uiStatusBar.showMessage(msg)

            if self.catalogue.type == "local":
                qimport.parse(self.catalogue.local_input, self.catalogue,
                              self.network)
            else:
                self.catalogue.sources_from_datacentre()
                self.catalogue.arrivals_from_datacentre(phases=["P", "S", "SKS"])

            if len(self.catalogue) == 0:
                msg = "No available sources from file/datacentre."
                self.uiStatusBar.showMessage(msg)

        self.uiSourceLoadButton.setEnabled(False)

        self._populateSourceList()

    def _populateCatalogueInformation(self, catalogue):
        """
        Populates the Catalogue Info widget with catalogue information

        Parameters
        ----------
        params : dict
            Dictionary object that contains the parameters that define the
            catalogue

        """

        self.uiStartDateDisplay.setText(catalogue.sdate)
        self.uiEndDateDisplay.setText(catalogue.edate)
        self.uiMinRadDisplay.setText(catalogue.minrad)
        self.uiMaxRadDisplay.setText(catalogue.maxrad)
        self.uiMinMagDisplay.setText(catalogue.minmag)
        self.uiSourceCountDisplay.setText(str(len(catalogue)))

    def plotCatalogueMap(self, replot=False):
        self.uiStatusBar.showMessage("Plotting catalogue map...")

        self.map_bg = None
        self.map_drag = False
        self.map_click = None

        self.uiMapMpl.show()

        if replot:
            self.uiMapMpl.canvas.ax.clear()
            lims = {"lon0": float(self.uiMinLonInput.text()),
                    "lon1": float(self.uiMaxLonInput.text()),
                    "lat0": float(self.uiMinLatInput.text()),
                    "lat1": float(self.uiMaxLatInput.text())}
        else:
            lims = None

        self.catalogue.plot(map_widget=self.uiMapMpl,
                            lims=lims)

        if not replot:
            if self.catalogue.type == "local":
                self.uiMinLonInput.setText(str(f"{self.catalogue.lon0:.5f}"))
                self.uiMaxLonInput.setText(str(f"{self.catalogue.lon1:.5f}"))
                self.uiMinLatInput.setText(str(f"{self.catalogue.lat0:.5f}"))
                self.uiMaxLatInput.setText(str(f"{self.catalogue.lat1:.5f}"))

        # Connect to the map to grab background once Qt has done resizing
        self.uiMapMpl.canvas.mpl_connect("draw_event", self._mapDrawEvent)

        # Draw the canvas once Qt has done all resizing
        self.uiMapMpl.canvas.draw_idle()

        self.uiStatusBar.showMessage("Catalogue map plot complete.")

    def _mapDrawEvent(self, event):
        """
        Grab the map background when it is drawn

        Parameters
        ----------
        event : DrawEvent
            Contains information about the draw event

        """

        self.map_bg = event.canvas.copy_from_bbox(event.canvas.ax.bbox)

    def _onPick(self, event):
        """
        Process picks on data in plot

        Parameters
        ----------
        event : PickEvent
            Contains information about the pick event

        """

        if event.mouseevent.key == "shift":
            return

        artist = event.artist

        label = artist.get_label()
        if "RECEIVER" in label:
            receiverid = label.split(": ")[1]
            receiver = self.catalogue.network.lookup(receiverid)
            self._updateReceiverInformation(receiver)

        if "SOURCE" in label:
            sourceid = label.split(": ")[1]
            source = self.catalogue.lookup(sourceid)
            self._updateSourceInformation(source)

            if event.mouseevent.key == "control":
                source.pick_window(filt=self.filt)

    def _onMapClick(self, event):
        """
        Grab and store the x and y positions of the click point

        """

        if event.key != "shift":
            return
        if self.catalogue.type == "teleseismic":
            return
        if self.map_drag:
            return
        if event.inaxes != self.uiMapMpl.canvas.ax:
            return

        xpress, ypress = event.xdata, event.ydata
        self.map_drag = True
        self.map_click = xpress, ypress

    def _onMapMove(self, event):
        """
        Draws a rectangle around map region to zoom if map_drag is True

        """

        if self.catalogue.type == "teleseismic":
            return
        if not self.map_drag:
            return
        if event.inaxes != self.uiMapMpl.canvas.ax:
            return

        # Update current position
        xpress, ypress = self.map_click
        x, y = event.xdata, event.ydata
        dx = x - xpress
        dy = y - ypress

        # Create Rectangle patch
        rect = Rectangle((xpress, ypress), dx, dy, edgecolor="red", fill=False,
                         transform=self.uiMapMpl.canvas.proj)

        # Add rectangle to map
        self.uiMapMpl.canvas.restore_region(self.map_bg)
        self.uiMapMpl.canvas.ax.add_patch(rect)
        self.uiMapMpl.canvas.ax.draw_artist(rect)
        self.uiMapMpl.canvas.blit(self.uiMapMpl.canvas.ax.bbox)

    def _onMapRelease(self, event):
        """
        Draws final rectangle around map region to zoom if map_drag if True

        """

        if self.catalogue.type == "teleseismic":
            return
        if not self.map_drag:
            return

        # Reuse code for movement
        self._onMapMove(event)

        # Update final position
        xpress, ypress = self.map_click
        x, y = event.xdata, event.ydata

        # Set the text values of the lon/lat input boxes
        self.uiMinLonInput.setText(str(f"{min(xpress, x):.2f}"))
        self.uiMaxLonInput.setText(str(f"{max(xpress, x):.2f}"))
        self.uiMinLatInput.setText(str(f"{min(ypress, y):.2f}"))
        self.uiMaxLatInput.setText(str(f"{max(ypress, y):.2f}"))

        # Reset map click and lock variables
        self.map_click = None
        self.map_drag = False

    def sourceSelect(self, index):
        """
        Updates the displayed Source information

        Parameters
        ----------
        index : int
            Index of the item picked from the Sources List widget

        """

        source = self.catalogue.lookup(
            self.uiSourcesListView.model().data(index))

        self._updateSourceInformation(source)

    def export(self, format_):
        # Test for existence of catalogue
        if not self.catalogue:
            return
        if format_ == "NLLoc":
            qexport.export_nonlinloc(self.catalogue)

    def sourcePick(self, index):
        """
        Opens a picking window for a given source

        If the waveform files have not been downloaded, it attempts to retrieve
        them from the archive. If the file does not exist, the Arrival is
        removed from the Source

        Parameters
        ----------
        index : int
            Index of the item picked from the Sources List widget

        """

        source = self.catalogue.lookup(
            self.uiSourcesListView.model().data(index))

        picker = Picker(self.catalogue, self.archive, self.network,
                        source, self.filt)

        # if not picker.exec_():
        #     print("bork")
        #     return
        # else:
        #     print("beep")
        #     if self.uiDefaultSaveAction.isChecked():
        #         self.catalogue.save()

    def _populateSourceList(self):
        """
        Populates the SourceList object with Source UIDs

        """

        sources = self.catalogue.sources

        model = QtGui.QStandardItemModel(self.uiSourcesListView)
        for uid in sources.keys():
            item = QtGui.QStandardItem(uid)
            model.appendRow(item)
        self.uiSourcesListView.setModel(model)

    def _updateReceiverInformation(self, rec):
        """
        Populates the Receiver Info widget with receiver information

        Parameters
        ----------
        rec : Quake Receiver object
            Describes the receiver

        """

        ddate = rec.deployment.isoformat().split("T")[0]
        rdate = rec.retrieval.isoformat().split("T")[0]

        self.uiReceiverNameDisplay.setText(rec.uid)
        self.uiReceiverLonDisplay.setText(f'{rec.longitude:.4f}')
        self.uiReceiverLatDisplay.setText(f'{rec.latitude:.4f}')
        self.uiReceiverElevDisplay.setText(f'{rec.elevation:.4f}')
        self.uiReceiverDepDisplay.setText(ddate)
        self.uiReceiverRetDisplay.setText(rdate)

    def _updateSourceInformation(self, src):
        """
        Populates the Source Info widget with source information

        Parameters
        ----------
        src : Quake Source object
            Describes the source

        """

        odate = src.otime.isoformat().split("T")[0]
        otime = src.otime.isoformat().split("T")[1].split(".")[0]

        self.uiOriginDateDisplay.setText(odate)
        self.uiOriginTimeDisplay.setText(otime)
        self.uiSourceLonDisplay.setText(f'{src.longitude:.4f}')
        self.uiSourceLatDisplay.setText(f'{src.latitude:.4f}')
        self.uiSourceDepthDisplay.setText(f'{src.depth:.4f}')
        if (type(src.magnitude) == float) or \
           (type(src.magnitude) == np.float64):
            self.uiSourceMagDisplay.setText(f'{src.magnitude:.2f}')
        else:
            self.uiSourceMagDisplay.setText(str(src.magnitude))


class NewCatalogueDialogue(qt.QDialog):
    """
    NewCatalogueDialogue class

    Dialogue window for generating a new catalogue

    Attributes
    ----------

    Methods
    -------

    """

    def __init__(self):
        """
        Class initialisation

        """

        super().__init__()

        self.path = None
        self.params = None
        self.archive = None
        self.network = None

        self.initUI()

    def initUI(self):
        """
        Load UI from designer file and initialise interface

        """

        uic.loadUi("gui/ui/new_catalogue.ui", self)

        self.connect()

        self.setWindowTitle("Quake - new catalogue")
        self.setWindowIcon(QtGui.QIcon("gui/ui/resources/icon.png"))
        self.show()

    def connect(self):
        """
        Connect all signals from Qt to functions

        """

        self.uiCatalogueDirectoryButton.clicked.connect(self.browseCatalogue)
        self.uiArchiveDirectoryButton.clicked.connect(self.browseArchive)
        self.uiReceiverFileButton.clicked.connect(self.browseReceivers)
        self.uiTeleseismicRadio.toggled.connect(self.teleRadioClicked)
        self.uiLocalRadio.toggled.connect(self.localRadioClicked)
        self.uiButtonBox.accepted.connect(self.actionAccept)
        self.uiButtonBox.rejected.connect(self.actionReject)
        self.uiLocalFileButton.clicked.connect(self.browseLocalFile)
        self.uiLocalPathButton.clicked.connect(self.browseLocalPath)
        self.uiLocalInputFormatComboBox.currentIndexChanged.connect(
            self.inputSelect)

    def browseCatalogue(self):
        """
        Opens a built in Qt window to select a catalogue directory

        """

        msg = "Choose catalogue directory"
        catalogue_path = qt.QFileDialog.getExistingDirectory(self, msg)
        self.uiCatalogueDirectoryInput.setText(catalogue_path)

    def browseArchive(self):
        """
        Opens a built in Qt window to select an archive directory

        """

        msg = "Choose archive directory"
        archive_path = qt.QFileDialog.getExistingDirectory(self, msg)
        self.uiArchiveDirectoryInput.setText(archive_path)

    def browseReceivers(self):
        """
        Opens a built in Qt window to select a receiver file

        """

        msg = "Open file"
        filename = qt.QFileDialog.getOpenFileName(self, msg)
        self.uiReceiverFileInput.setText(filename[0])

    def browseLocalFile(self):
        """
        Opens a built in Qt window to select a file containing earthquake
        sources

        """

        msg = "Open file"
        filename = qt.QFileDialog.getOpenFileName(self, msg)
        self.uiLocalFileInput.setText(filename[0])

    def browseLocalPath(self):
        """
        Opens a built in Qt window to select a directory containing earthquake
        sources

        """

        msg = "Choose QuakeMigrate directory"
        pathname = qt.QFileDialog.getExistingDirectory(self, msg)
        self.uiLocalPathInput.setText(pathname)

    def teleRadioClicked(self, enabled):
        """Control for teleseismic toggle"""

        if enabled:
            self.uiInputTypeStacked.setCurrentIndex(0)

    def localRadioClicked(self, enabled):
        """Control for local toggle"""

        if enabled:
            self.uiInputTypeStacked.setCurrentIndex(1)

    def inputSelect(self):
        """Control for handling different input formats"""

        if self.uiLocalInputFormatComboBox.currentText() == "QuakeMigrate":
            self.uiLocalOptionsStacked.setCurrentIndex(1)
        else:
            self.uiLocalOptionsStacked.setCurrentIndex(0)

    def actionAccept(self):
        """
        Accept the Dialogue

        Performs a series of checks to ensure the input form has been completed
        correctly

        """

        if (self.uiCatalogueDirectoryInput.text() == "") or \
           (self.uiCatalogueNameInput.text() == "") or \
           (self.uiReceiverFileInput.text() == ""):
            msg = "Please complete the form and try again."
            qt.QMessageBox.about(self, "Error!", msg)
            return

        # --- Path ---
        path = pathlib.Path(self.uiCatalogueDirectoryInput.text())
        if not path.is_dir():
            msg = "Invalid path to catalogue directory."
            qt.QMessageBox.about(self, "Error!", msg)
            return
        uid = self.uiCatalogueNameInput.text()

        self.path = "{}/{}".format(path, uid)

        if self.uiTeleseismicRadio.isChecked():
            type_ = "teleseismic"
            source = "IRIS"
            d = {"minmag": self.uiMinMagInput.text(),
                 "minrad": self.uiMinRadInput.text(),
                 "maxrad": self.uiMaxRadInput.text(),
                 "sdate": self.uiStartDateInput.date().toString(Qt.ISODate),
                 "edate": self.uiEndDateInput.date().toString(Qt.ISODate)}
        elif self.uiLocalRadio.isChecked():
            type_ = "local"
            source = self.uiLocalInputFormatComboBox.currentText()
            if source == "QuakeMigrate" or source == ".nonlinloc":
                local_input = pathlib.Path(self.uiLocalPathInput.text())
                if not local_input.is_dir():
                    msg = "Invalid path to input directory."
                    qt.QMessageBox.about(self, "Error!", msg)
                    return
            else:
                local_input = pathlib.Path(self.uiLocalFileInput.text())
                if not local_input.is_file():
                    msg = "Invalid path to input file."
                    qt.QMessageBox.about(self, "Error!", msg)
                    return
            d = {"local_input": str(local_input)}
        params = {"uid": uid,
                  "type": type_,
                  "source": source,
                  "cdate": QDateTime.currentDateTime().toUTC().toString(Qt.ISODate)}

        self.params = {**params, **d}

        # --- Archive ---
        archive = pathlib.Path(self.uiArchiveDirectoryInput.text())
        if not archive.is_dir():
            msg = "Invalid path to archive directory."
            qt.QMessageBox.about(self, "Error!", msg)
            return
        archive_schema = self.uiArchiveFormatComboBox.currentText()
        if archive_schema == "custom":
            structure = self.uiArchiveFormatInput.text()
        else:
            structure = ""

        self.archive = {"path": str(archive),
                        "schema": archive_schema,
                        "structure": structure}

        # --- Network ---
        receivers = pathlib.Path(self.uiReceiverFileInput.text())
        if not receivers.is_file():
            msg = "Invalid path to receiver file."
            qt.QMessageBox.about(self, "Error!", msg)
            return

        self.network = str(receivers)

        self.accept()

    def actionReject(self):
        """Send reject signal to Dialog"""

        self.reject()


class TelePhaseDialogue(qt.QDialog):

    def __init__(self):
        super().__init__()

        self.phases = []

        self.initUI()

    def initUI(self):
        uic.loadUi("gui/ui/teleseismic_phases.ui", self)

        self.connect()

        # Populate the receiver list
        self.model = QtGui.QStandardItemModel(self.uiPhaseListView)

        self.setWindowTitle("Quake - teleseismic phase selection")
        self.setWindowIcon(QtGui.QIcon("gui/ui/resources/icon.png"))
        self.show()

    def connect(self):
        self.uiAddPhaseButton.clicked.connect(self.addPhase)
        self.uiButtonBox.accepted.connect(self.actionAccept)
        self.uiButtonBox.rejected.connect(self.actionReject)

    def addPhase(self):
        # Grab the current phase to be added
        phase = self.uiInputPhaseComboBox.currentText()

        # Check if phase already in list
        if phase in self.phases:
            pass
        # If not, add it
        else:
            self.phases.append(phase)

            # Add it to the ListView object
            item = QtGui.QStandardItem(phase)
            self.model.appendRow(item)

            # Update model
            self.uiPhaseListView.setModel(self.model)

    def actionAccept(self):
        if not self.phases:
            qt.QMessageBox.about(self, "Error!", "No phases selected.")
            return

        # Update the input phases in the parent class
        self.parent.input_phases = self.phases

        # Send accept signal to Dialog
        self.accept()

    def actionReject(self):
        # Send reject signal to Dialog
        self.reject()


class DefaultFilterDialogue(qt.QDialog):

    def __init__(self, parent):
        super().__init__()

        self.parent = parent

        self.initUI()

    def initUI(self):
        uic.loadUi("gui/ui/default_filter.ui", self)

        # Connect to actions and interactive widgets
        self.connect()

        self.setWindowTitle("Quake - default filter")
        self.setWindowIcon(QtGui.QIcon("gui/ui/resources/icon.png"))
        self.show()

    # ------------------------------

    # ------------------------------
    # Connection and event functions
    # ------------------------------

    def connect(self):
        # Connect to accept/reject actions
        self.uiButtonBox.accepted.connect(self.actionAccept)
        self.uiButtonBox.rejected.connect(self.actionReject)

    # ------------------------------

    # ------------------------------
    # Accept/Reject action overrides
    # ------------------------------

    def actionAccept(self):
        # Builds a dictionary of keyword arguments for the filter
        # function in the Source class, based on filter type
        filter_type = self.uiFilterTypeComboBox.currentText()
        corners = int(self.uiPoleCountComboBox.currentText())
        zerophase = self.uiZeroPhaseCheckBox.isChecked()

        filter_options = {"type": filter_type,
                          "corners": corners,
                          "zerophase": zerophase}

        if filter_type == "bandpass":
            try:
                freqmin = float(self.uiLowFreqInput.text())
            except ValueError:
                qt.QMessageBox.about(self, "Error!", "Invalid lowpass frequency.")
                return

            try:
                freqmax = float(self.uiHighFreqInput.text())
            except ValueError:
                qt.QMessageBox.about(self, "Error!", "Invalid highpass frequency.")
                return

            # Create filter options dictionary
            filter_freqs = {"freqmin": freqmin,
                            "freqmax": freqmax}

        elif filter_type == "lowpass":
            try:
                freq = float(self.uiLowFreqInput.text())
            except ValueError:
                qt.QMessageBox.about(self, "Error!", "Invalid lowpass frequency.")
                return

            # Create filter options dictionary
            filter_freqs = {"freq": freq}

        elif filter_type == "highpass":
            try:
                freq = float(self.uiHighFreqInput.text())
            except ValueError:
                qt.QMessageBox.about(self, "Error!", "Invalid highpass frequency.")
                return

            # Create filter options dictionary
            filter_freqs = {"freq": freq}

        # Update the default filter in parent class
        self.parent.filt = {**filter_options, **filter_freqs}

        # Send accept signal to Dialog
        self.accept()

    def actionReject(self):
        # Send reject signal to Dialog
        self.reject()


class WadatiWindow(qt.QMainWindow):

    wad_labels = {"x": "P traveltime / s",
                  "y": "S - P traveltime / s"}

    dist_labels = {"x": "Epicentral distance / km",
                   "y": "Traveltime / s"}

    def __init__(self, parent):
        super().__init__()

        self.ptravels = []
        self.stravels = []
        self.sptimes = []
        self.receivers = []
        self.dists = []

        self.parent = parent

        self.initUI()

    def initUI(self):
        uic.loadUi("gui/ui/wadati_window.ui", self)

        self.uiWadatiPlotMpl.canvas.tracePlot(labels=self.wad_labels)
        self.uiDistancePlotMpl.canvas.tracePlot(labels=self.dist_labels)

        self.connect()

        self.setWindowTitle("Quake - Wadati plot window")
        self.setWindowIcon(QtGui.QIcon("gui/ui/resources/icon.png"))
        self.show()

    def connect(self):
        # Calculate Vp/Vs button
        self.uiCalcVpVsButton.clicked.connect(self.calcVpVs)

        # Plot connection
        self.uiWadatiPlotMpl.canvas.mpl_connect("pick_event", self._onPick)
        self.uiDistancePlotMpl.canvas.mpl_connect("pick_event", self._onPick)

    def calcVpVs(self):
        # Insert code to calc best fitting line to p travel and s-p times
        # Should be a straight line
        pass

    def _onPick(self, event):
        artist = event.artist
        label = artist.get_label()
        self.uiReceiverDisplay.setText(label)

    def addPick(self, ptravel, stravel, dist, receiver):
        # Calculate sptime
        sptime = stravel - ptravel
        # Add a p traveltime and SP-time to the plot
        self.ptravels.append(ptravel)
        self.stravels.append(stravel)
        self.sptimes.append(sptime)
        self.dists.append(dist)
        self.receivers.append(receiver)

        self.plotWadati(receiver)
        self.plotDistance(receiver)

    def plotWadati(self, receiver):
        wadati_canvas = self.uiWadatiPlotMpl.canvas

        tolerance = 10
        for i in range(len(self.ptravels)):
            wadati_canvas.ax.scatter(self.ptravels[i], self.sptimes[i], 12,
                                     marker="o", color="k", picker=tolerance,
                                     zorder=10, label=self.receivers[i])

        wadati_canvas.draw_idle()

    def plotDistance(self, receiver):
        distance_canvas = self.uiDistancePlotMpl.canvas

        tolerance = 10
        for i in range(len(self.ptravels)):
            distance_canvas.ax.scatter(self.dists[i], self.ptravels[i], 12,
                                       marker="o", color="red", picker=tolerance,
                                       zorder=10, label=self.receivers[i])
            distance_canvas.ax.scatter(self.dists[i], self.stravels[i], 12,
                                       marker="o", color="blue", picker=tolerance,
                                       zorder=10, label=self.receivers[i])

        distance_canvas.draw_idle()

    # Function from TG - uses orthogonal distance regression
    def fitLine(self, p_times, sp_times, p_errors=None, s_errors=None):
        def func(p, x):
            x = np.array(x)
            return p[0] * x + p[1]
        linear = odr.Model(func)
        my_data = odr.RealData(p_times, sp_times, sx=p_errors, sy=sp_errors)
        my_odr  = odr.ODR(my_data, linear, beta0=[0.75, 0.])
        results = my_odr.run()


if __name__ == "__main__":

    app = qt.QApplication(sys.argv)

    # Set stylesheet
    file = QFile("gui/ui/style.qss")
    file.open(QFile.ReadOnly | QFile.Text)
    stream = QTextStream(file)
    app.setStyleSheet(stream.readAll())

    quake = Quake()
    sys.exit(app.exec_())
