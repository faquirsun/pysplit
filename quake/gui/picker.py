# -*- coding: utf-8 -*-
"""
This module requires that "numpy", "pandas", "matplotlib" and "PyQt5" are
installed in your Python environment

The back-end is powered by a series of classes in quake.core

Author: Conor Bacon
"""

from matplotlib.patches import Rectangle
import numpy as np
from PyQt5.QtCore import Qt
import PyQt5.QtWidgets as qt
import PyQt5.QtGui as QtGui
from PyQt5 import uic

from quake.core.pick import Pick


class Picker(qt.QMainWindow):
    """
    Picker class

    Provides a GUI for picking the phase arrivals of seismic waves.

    """

    def __init__(self, catalogue, archive, network, source, filt=None):
        """
        Class initialisation method

        Parameters
        ----------
        source : Quake Source object
            Contains information about a seismic source

        """

        super().__init__()

        self.catalogue = catalogue
        self.archive = archive
        self.network = network
        self.source = source
        self.filt = filt
        self.default_filt = filt

        # Initialise attribute defaults
        self.pick_type = "P"
        self.pick_line_color = "red"
        self.lines = {}
        self.components = "ZNE"
        self.show_spectrogram = False

        self.arrivals = iter(source)

        try:
            uid, self.arrival = next(self.arrivals)
            self.arrival.load_waveform(self.archive, self.catalogue.path)
            self.initUI()
        except StopIteration:
            print("No available arrivals.")
            self.close()

    def initUI(self):
        """
        Load UI from designer file and initialise interface

        """

        uic.loadUi("gui/ui/pick.ui", self)

        self.showMaximized()

        self.connect()

        self.uiComponentsMpl.canvas.tracePlot()
        self._updateSourceInformation(self.arrival.source)
        self._updateReceiverInformation(self.arrival.receiver)
        if self.default_filt is not None:
            self._populateFilter()
        self.plotTraces()

        self.setWindowTitle("Quake - trace picker")
        self.setWindowIcon(QtGui.QIcon("gui/ui/resources/icon.png"))
        self.show()

    def connect(self):
        """
        Connect all signals from Qt to functions

        """

        self.uiSaveAction.triggered.connect(self.save)

        self.uiPPickAction.triggered.connect(lambda: self.updatePick("P"))
        self.uiSPickAction.triggered.connect(lambda: self.updatePick("S"))
        self.uiCustomPickAction.triggered.connect(lambda: self.updatePick("C"))

        self.uiApplyFilterButton.clicked.connect(self.applyFilter)
        self.uiRemoveFilterButton.clicked.connect(self.removeFilter)

        self.uiRotateCompButton.clicked.connect(self.rotateComponents)
        # self.uiViewSpectrogramButton.clicked.connect(self.toggleSpectrogram)
        self.uiNextTraceButton.clicked.connect(
            lambda: self.changeTrace(move=1))
        self.uiLastTraceButton.clicked.connect(
            lambda: self.changeTrace(move=-1))
        self.uiRejectTraceButton.clicked.connect(self.rejectTrace)
        self.uiToggleLimitsButton.clicked.connect(self.toggleLims)
        self.uiResetPlotButton.clicked.connect(self.resetPlot)

    def save(self):
        """
        Save the current state of the catalogue, including any picks made
        since the last save

        """

        self.catalogue.save()

    def _plotconnect(self):
        """
        Connect to plot callbacks.

        """

        canvas = self.uiComponentsMpl.canvas

        self._mid = canvas.mpl_connect("motion_notify_event", self._onMove)
        self._cid = canvas.mpl_connect("button_press_event", self._onClick)
        self._rid = canvas.mpl_connect("button_release_event", self._onRelease)

    def _plotdisconnect(self):
        """
        Temporarily disconnect current plot callbacks to avoid recursion.

        """

        canvas = self.uiComponentsMpl.canvas

        canvas.mpl_disconnect(self._mid)
        canvas.mpl_disconnect(self._cid)
        canvas.mpl_disconnect(self._rid)

    def keyPressEvent(self, event):
        """
        Override default callback for KeyPressEvent

        Parameters
        ----------
        event : KeyPressEvent
            Contains information about the key press event

        """

        if event.key() == Qt.Key_U or event.key() == Qt.Key_D:
            if self.pick_type not in self.lines.keys():
                return
            else:
                pick_polarity = chr(event.key())
                self.uiPolarityDisplay.setText(pick_polarity)
                self.arrival.picks[self.pick_type].polarity = pick_polarity

    def _onMove(self, event):
        """
        Callback for motion notify event

        Parameters
        ----------
        event : motion event
            Contains information about the mouse motion event

        """

        canvas = self.uiComponentsMpl.canvas

        if event.inaxes not in canvas.axes:
            return
        if event.key == "shift":
            if not self.trace_drag:
                return

            for bg in self.bgs:
                canvas.restore_region(bg)

            xpress, ypress = self.trace_click
            xmove, ymove = event.xdata, event.ydata
            dx = xmove - xpress
            dy = ymove - ypress

            for ax in canvas.axes:
                rect = Rectangle((xpress, ypress), dx, dy,
                                 edgecolor="red", fill=False)
                ax.add_patch(rect)
                ax.draw_artist(rect)

        else:
            for bg in self.bgs:
                canvas.restore_region(bg)

            x = event.xdata

            for ax, cursor in zip(canvas.axes, self.cursors):
                cursor.set_data([x, x], self.y)
                ax.draw_artist(cursor)

        self._plotPicks()

        for ax in canvas.axes:
            canvas.blit()

    def _onClick(self, event):
        """
        Callback for mouse click event

        Parameters
        ----------
        event : click event
            Contains information about the mouse click event

        """

        if event.key == "shift":
            if self.trace_drag:
                return
            self.trace_drag = True
            self.trace_click = event.xdata, event.ydata
            return

        canvas = self.uiComponentsMpl.canvas

        adjusted_xdata = round(event.xdata / self.arrival.delta) \
            * self.arrival.delta

        if event.button == 1:
            color = self.pick_line_color
            pick_type = self.pick_type

            pick = Pick(pick_type, adjusted_xdata
                        - (self.source.otime - self.arrival.starttime))
            self.arrival + pick

            pick_time = self.arrival.starttime + adjusted_xdata
            self.uiPickTimeDisplay.setText(pick_time.isoformat())
            self.uiPhaseDisplay.setText(pick_type)
        elif event.button == 3:
            # Implement pick delete option
            return

        lines = []
        for ax in canvas.axes:
            lines.append(ax.axvline(adjusted_xdata, linewidth=1,
                                    color=color, animated=True))

        self.lines[pick_type] = lines

        for bg in self.bgs:
            canvas.restore_region(bg)

        self._plotPicks()

        for ax in canvas.axes:
            canvas.blit(ax.bbox)

    def _onRelease(self, event):
        canvas = self.uiComponentsMpl.canvas

        if event.inaxes not in canvas.axes:
            return
        if event.key != "shift":
            return
        if not self.trace_drag:
            return

        xpress, ypress = self.trace_click
        xrelease, yrelease = event.xdata, event.ydata

        xlims = (min(xpress, xrelease), max(xpress, xrelease))
        ylims = (min(ypress, yrelease), max(ypress, yrelease))

        self.previous_lims = self.lims
        self.lims = [xlims, ylims]

        self.plotTraces(replot=True)

    def plotTraces(self, replot=False):
        canvas = self.uiComponentsMpl.canvas

        for ax in canvas.axes:
            ax.clear()

        # Create background variable
        self.bgs = []
        self.cursors = []

        self.trace_drag = False
        self.trace_click = None
        try:
            if not replot:
                # Add limits variable
                self.previous_lims = None
                self.lims = None

                # Set filter to default filter
                self.filt = self.default_filt

                # Look up any picks and make vertical lines
                for uid, pick in self.arrival.picks.items():
                    pick_times = pick.ttime
                    c = "gray"
                    for time in pick_times:
                        lines = []
                        for ax in canvas.axes:
                            line = ax.axvline(time
                                              + (self.source.otime
                                                 - self.arrival.starttime),
                                              linewidth=0.75,
                                              color=c, alpha=0.4)
                            lines.append(line)

                            self.lines[uid] = lines

            self.arrival.filter(filt=self.filt)

            # Plot the traces
            self.arrival.plot(canvas.axes, lims=self.lims)
        except AttributeError:
            pass

        # Connect to trace to grab the background once Qt has done resizing
        canvas.mpl_connect("draw_event", self._drawEvent)

        # Draw when Qt has done all resizing
        canvas.draw_idle()

        # Initialise the cursor to track mouse position on the axes
        for ax in canvas.axes:
            self.cursors.append(ax.axvline(5, linewidth=1, color="0.5",
                                           animated=True))
        _, self.y = self.cursors[0].get_data()

        self._plotconnect()

    def _drawEvent(self, event):
        """
        Grab the trace background when it is drawn

        Parameters
        ----------
        event : drawEvent
            Contains information about the draw event

        """

        for ax in event.canvas.axes:
            self.bgs.append(event.canvas.copy_from_bbox(ax.bbox))

    def _plotPicks(self):
        """
        Replot any picked lines

        """

        canvas = self.uiComponentsMpl.canvas

        for pick_type, lines in self.lines.items():
            for ax in canvas.axes:
                for line in lines:
                    ax.draw_artist(line)

    def _updateReceiverInformation(self, rec):
        """
        Populates the Receiver Info widget with receiver information

        Parameters
        ----------
        rec : Quake Receiver object
            Describes the receiver

        """

        self.uiReceiverNameDisplay.setText(rec.uid)
        self.uiReceiverLonDisplay.setText(f'{rec.longitude:.4f}')
        self.uiReceiverLatDisplay.setText(f'{rec.latitude:.4f}')
        self.uiReceiverElevDisplay.setText(f'{rec.elevation:.4f}')

        try:
            ddate = rec.deployment.isoformat().split("T")[0]
            rdate = rec.retrieval.isoformat().split("T")[0]

            self.uiReceiverDepDisplay.setText(ddate)
            self.uiReceiverRetDisplay.setText(rdate)
        except AttributeError:
            self.uiReceiverDepDisplay.setText("n/a")
            self.uiReceiverRetDisplay.setText("n/a")

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

    def applyFilter(self):
        """
        Filter the seismic traces

        """

        try:
            filt_type = self.uiFilterTypeComboBox.currentText()
            corners = int(self.uiPoleCountComboBox.currentText())
            zerophase = self.uiZeroPhaseCheckBox.isChecked()
            freqmin = float(self.uiLowFreqInput.text())
            freqmax = float(self.uiHighFreqInput.text())

            self.filt = {"type": filt_type,
                         "corners": corners,
                         "zerophase": zerophase,
                         "freqmin": freqmin,
                         "freqmax": freqmax}
        except ValueError:
            msg = "Fill in all filter options correctly."
            qt.QMessageBox.about(self, "Warning!", msg)
            return

        self._plotdisconnect()
        self.plotTraces(replot=True)

    def removeFilter(self):
        """
        Remove any filters applied to the current Arrival

        """

        self._plotdisconnect()
        self.filt = None
        self.arrival.filter()
        self.plotTraces(replot=True)

    def _populateFilter(self):
        """
        Populate filter input with filter parameters

        """

        filter_type = self.default_filt["type"]

        self.uiFilterTypeComboBox.setCurrentText(filter_type)
        self.uiPoleCountComboBox.setCurrentText(str(self.default_filt["corners"]))

        if self.default_filt["zerophase"]:
            self.uiZeroPhaseCheckBox.setChecked(True)

        if filter_type == "bandpass":
            self.uiLowFreqInput.setText(str(self.default_filt["freqmin"]))
            self.uiHighFreqInput.setText(str(self.default_filt["freqmax"]))

        if filter_type == "lowpass":
            self.uiLowFreqInput.setText(str(self.default_filt["freq"]))

        if filter_type == "highpass":
            self.uiHighFreqInput.setText(str(self.default_filt["freq"]))

    def updatePick(self, pick_type):
        """
        Updates pick type and corresponding variables

        Parameters
        ----------
        pick_type : str
            String specifying the type of pick ("S", "P", "C")

        """

        if pick_type == "P":
            self.pick_type = pick_type
            self.pick_line_color = "red"

        elif pick_type == "S":
            self.pick_type = pick_type
            self.pick_line_color = "blue"

        elif pick_type == "C":
            # Only run if the radio button is being toggled on
            if not self.uiCustomRadio.isChecked():
                return
            dlg = CustomPickDialogue()

            if not dlg.exec_():
                return

            self.pick_type = dlg.pick_type
            self.pick_line_color = dlg.pick_line_color

            self.uiCustomPhaseDisplay.setText(self.pick_type)

        self.uiPolarityDisplay.setText("")

    def rotateComponents(self):
        """
        Rotate Arrival components

        """

        if self.components == "ZNE":
            self.components = "ZRT"
            self.arrival.rotate(method="NE->RT")

        elif self.components == "ZRT":
            self.components = "ZNE"
            self.arrival.rotate(method="RT->NE")

        self.plotTraces(replot=True)

    def changeTrace(self, move):
        """
        Moves the trace counter and updates the current source-receiver pair.

        Parameters
        ----------
        move : int
            Value to change the counter by (-1, 0 or 1)

        """

        # self.counter += move

        # if self.counter == -1:
        #     self.counter = len(self.sr_pairs) - 1
        # elif self.counter == len(self.sr_pairs):
        #     self.counter = 0

        try:
            uid, self.arrival = next(self.arrivals)
        except StopIteration:
            self.arrivals = iter(self.source)
            uid, self.arrival = next(self.arrivals)

        self.arrival.load_waveform(self.archive, self.catalogue.path)

        self._updateSourceInformation(self.arrival.source)
        self._updateReceiverInformation(self.arrival.receiver)
        self._updateTrace()

    def _updateTrace(self):
        # Check if P and S have been picked (and window exists)
        if "P_manual" in self.arrival.picks.keys() and \
           "S_manual" in self.arrival.picks.keys() and \
           self.wadatiWindow:
            ptravel = self.arrival.starttime - self.arrival.otime
            ptravel += self.arrival.picks["P_manual"]["ttime"]
            stravel = self.arrival.starttime - self.arrival.otime
            stravel += self.arrival.picks["S_manual"]["ttime"]
            self.wadatiWindow.addPick(ptravel, stravel, self.receiver)

        self.uiPRadio.setChecked(True)

        self.lines = {}
        self.pick_time = ""
        self.pick_phase = ""

        self.uiPolarityDisplay.setText("")
        self.uiPickTimeDisplay.setText(self.pick_time)
        self.uiPhaseDisplay.setText(self.pick_phase)

        try:
            self._plotdisconnect()
        except AttributeError:
            pass
        self.plotTraces()

    def rejectTrace(self):
        """
        Removes the current source-receiver pair from the arrivals catalogue.

        Unlinks the source-receiver MSEED files and deletes the corresponding
        entry in

        """

        print("Deleting trace...")
        self.source - self.arrival

        del self.sr_pairs[self.counter]

        self.changeTrace(move=0)

    def toggleLims(self):
        """
        Temporary function to implement a simple undo/redo for the plot
        details.

        TO-DO
        -----
        Implement a proper system of undo/redo

        """

        self.lims, self.previous_lims = self.previous_lims, self.lims

        self.plotTraces(replot=True)

    def resetPlot(self):
        """
        Reset the plot of the traces, returning to the default limits and
        filter.

        """

        self._plotdisconnect()
        self.lims = None
        self.filt = self.default_filt
        self.plotTraces()


class CustomPickDialogue(qt.QDialog):
    """
    CustomPickDialogue class

    Dialogue window for selecting a custom pick type

    """

    def __init__(self):
        """
        Class initialisation

        """

        super().__init__()

        self.initUI()

    def initUI(self):
        """
        Load UI from designer file and initialise interface

        """

        uic.loadUi("gui/ui/custom_phase.ui", self)

        self.setWindowTitle("Quake - custom phase pick")
        self.setWindowIcon(QtGui.QIcon("gui/ui/resources/icon.png"))
        self.show()

    def actionAccept(self):
        """
        Accept the Dialogue

        Performs a series of checks to ensure the input form has been completed
        correctly

        """

        try:
            self.pick_type = self.uiPhaseComboBox.currentText()
            self.pick_line_color = "orange"
        except ValueError:
            qt.QMessageBox.about(self, "Error!", "You need to specify a phase to pick!")

        self.accept()

    def actionReject(self):
        """Send reject signal to Dialogue"""

        self.reject()
