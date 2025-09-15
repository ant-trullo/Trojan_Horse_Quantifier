"""This is the main window of the software to analyze colocalization between VSV and delta viruses.
This is version 1.0, since September 2025.

author antonio.trullo@igmm.cnrs.fr
"""

from importlib import reload
import sys
import traceback
import numpy as np
import pyqtgraph as pg
from PyQt6 import QtGui, QtWidgets, QtCore
from PyQt6.QtGui import QAction
from PyQt6.QtWidgets import QVBoxLayout

import LoadRawData
import DeltaDetector
import VsvDetector
import AnalysisSaver


class MainWindow(QtWidgets.QMainWindow):
    """Main windows: coordinates all the actions, algorithms, visualization tools and analysis tools."""
    def __init__(self, parent=None):

        QtWidgets.QMainWindow.__init__(self, parent)

        ksf_h  =  np.load('keys_size_factor.npy')[0]
        ksf_w  =  np.load('keys_size_factor.npy')[1]

        widget  =  QtWidgets.QWidget(self)
        self.setCentralWidget(widget)

        load_data_action  =  QAction(QtGui.QIcon('Icons/load-hi.png'), "&Load Data", self)
        load_data_action.setShortcut("Ctrl+L")
        load_data_action.setStatusTip("Load raw data files")
        load_data_action.triggered.connect(self.load_raw_data)

        save_analysis_action  =  QAction(QtGui.QIcon('Icons/save-md.png'), "&Save Analysis", self)
        save_analysis_action.setShortcut("Ctrl+S")
        save_analysis_action.setStatusTip("Save Analysis")
        save_analysis_action.triggered.connect(self.save_analysis)

        load_analysis_action  =  QAction(QtGui.QIcon('Icons/save-md.png'), "&Load Analysis", self)
        load_analysis_action.setShortcut("Ctrl+W")
        load_analysis_action.setStatusTip("Load Analysis")
        load_analysis_action.triggered.connect(self.load_analysis)

        post_statistics_action  =  QAction(QtGui.QIcon('Icons/summarize.png'), "&Summarize Analyses", self)
        post_statistics_action.setShortcut("Ctrl+M")
        post_statistics_action.setStatusTip("Save Analysis")
        post_statistics_action.triggered.connect(self.post_statistics)

        exit_action  =  QAction(QtGui.QIcon('Icons/exit.png'), "&Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.setStatusTip("Exit application")
        exit_action.triggered.connect(self.close)

        menubar   =  self.menuBar()

        file_menu  =  menubar.addMenu("&File")
        file_menu.addAction(load_data_action)
        file_menu.addAction(save_analysis_action)
        file_menu.addAction(load_analysis_action)
        file_menu.addAction(exit_action)

        postprocessing_menu  =  menubar.addMenu("&PostProcessing")
        postprocessing_menu.addAction(post_statistics_action)

        # ~~~~~~~~~~~ VSV START ~~~~~~~~~~~~ #
        tabs_vsv  =  QtWidgets.QTabWidget()
        tab1_vsv  =  QtWidgets.QWidget()
        tab2_vsv  =  QtWidgets.QWidget()

        frame_vsv_raw  =  pg.ImageView(self, name="FrameVSVRaw")
        frame_vsv_raw.ui.roiBtn.hide()
        frame_vsv_raw.ui.menuBtn.hide()
        frame_vsv_raw.view.setXLink("FrameVSVSegm")
        frame_vsv_raw.view.setYLink("FrameVSVSegm")

        frame_vsv_segm  =  pg.ImageView(self, name="FrameVSVSegm")
        frame_vsv_segm.ui.roiBtn.hide()
        frame_vsv_segm.ui.menuBtn.hide()
        frame_vsv_segm.view.setXLink("FrameVSVRaw")
        frame_vsv_segm.view.setYLink("FrameVSVRaw")

        frame_vsv_raw_box  =  QtWidgets.QHBoxLayout()
        frame_vsv_raw_box.addWidget(frame_vsv_raw)

        frame_vsv_segm_box  =  QtWidgets.QHBoxLayout()
        frame_vsv_segm_box.addWidget(frame_vsv_segm)

        tab1_vsv.setLayout(frame_vsv_raw_box)
        tab2_vsv.setLayout(frame_vsv_segm_box)

        tabs_vsv.addTab(tab1_vsv, "Raw Data")
        tabs_vsv.addTab(tab2_vsv, "Segmented")

        fname_vsv_lbl  =  QtWidgets.QLabel("VSV   ", self)

        segment_vsv_btn  =  QtWidgets.QPushButton("Segment", self)
        segment_vsv_btn.setToolTip("Segment vsv virus")
        segment_vsv_btn.clicked.connect(self.segment_vsv)
        segment_vsv_btn.setFixedSize(int(ksf_h * 130), int(ksf_w * 25))
        segment_vsv_btn.setEnabled(True)

        filter_vsv_btn  =  QtWidgets.QPushButton("Filter", self)
        filter_vsv_btn.clicked.connect(self.filter_vsv)
        filter_vsv_btn.setToolTip("Filter segmented vsv viruses")
        filter_vsv_btn.setFixedSize(int(ksf_h * 130), int(ksf_w * 25))
        filter_vsv_btn.setEnabled(True)

        thickness_thr_lbl  =  QtWidgets.QLabel('Thick Thr', self)
        thickness_thr_lbl.setFixedSize(int(ksf_h * 80), int(ksf_w * 25))

        thickness_thr_edt  =  QtWidgets.QLineEdit(self)
        thickness_thr_edt.textChanged[str].connect(self.thickness_thr_var)
        thickness_thr_edt.setToolTip("Thickness Threshold (in pixels)")
        thickness_thr_edt.setFixedSize(int(ksf_h * 35), int(ksf_w * 25))
        thickness_thr_edt.setText("18.")

        thickness_thr_box  =  QtWidgets.QHBoxLayout()
        thickness_thr_box.addWidget(thickness_thr_lbl)
        thickness_thr_box.addWidget(thickness_thr_edt)

        up_ratio_thr_lbl  =  QtWidgets.QLabel('Up Ratio', self)
        up_ratio_thr_lbl.setFixedSize(int(ksf_h * 80), int(ksf_w * 25))

        up_ratio_thr_edt  =  QtWidgets.QLineEdit(self)
        up_ratio_thr_edt.textChanged[str].connect(self.up_ratio_thr_var)
        up_ratio_thr_edt.setToolTip("Upper ratio threshold between major and minor axis of the ellipsoidal fitting")
        up_ratio_thr_edt.setFixedSize(int(ksf_h * 35), int(ksf_w * 25))
        up_ratio_thr_edt.setText("5.")

        up_ratio_thr_box  =  QtWidgets.QHBoxLayout()
        up_ratio_thr_box.addWidget(up_ratio_thr_lbl)
        up_ratio_thr_box.addWidget(up_ratio_thr_edt)

        low_ratio_thr_lbl  =  QtWidgets.QLabel('Low Ratio', self)
        low_ratio_thr_lbl.setFixedSize(int(ksf_h * 80), int(ksf_w * 25))

        low_ratio_thr_edt  =  QtWidgets.QLineEdit(self)
        low_ratio_thr_edt.textChanged[str].connect(self.low_ratio_thr_var)
        low_ratio_thr_edt.setToolTip("Upper ratio threshold between major and minor axis of the ellipsoidal fitting")
        low_ratio_thr_edt.setFixedSize(int(ksf_h * 35), int(ksf_w * 25))
        low_ratio_thr_edt.setText("1.1")

        low_ratio_thr_box  =  QtWidgets.QHBoxLayout()
        low_ratio_thr_box.addWidget(low_ratio_thr_lbl)
        low_ratio_thr_box.addWidget(low_ratio_thr_edt)

        high_area_thr_lbl  =  QtWidgets.QLabel('Area Up Thr', self)
        high_area_thr_lbl.setFixedSize(int(ksf_h * 80), int(ksf_w * 25))

        high_area_thr_edt  =  QtWidgets.QLineEdit(self)
        high_area_thr_edt.textChanged[str].connect(self.high_area_thr_var)
        high_area_thr_edt.setToolTip("Remove big objects that cannot be viruses")
        high_area_thr_edt.setFixedSize(int(ksf_h * 35), int(ksf_w * 25))
        high_area_thr_edt.setText("150")

        high_area_thr_box  =  QtWidgets.QHBoxLayout()
        high_area_thr_box.addWidget(high_area_thr_lbl)
        high_area_thr_box.addWidget(high_area_thr_edt)

        keys_vsv  =  QtWidgets.QVBoxLayout()
        keys_vsv.addWidget(segment_vsv_btn)
        keys_vsv.addLayout(thickness_thr_box)
        keys_vsv.addLayout(up_ratio_thr_box)
        keys_vsv.addLayout(low_ratio_thr_box)
        keys_vsv.addLayout(high_area_thr_box)
        keys_vsv.addWidget(filter_vsv_btn)
        keys_vsv.addStretch()

        tab_fname_vsv_box  =  QVBoxLayout()
        tab_fname_vsv_box.addWidget(fname_vsv_lbl)
        tab_fname_vsv_box.addWidget(tabs_vsv)

        layout_vsv  =  QtWidgets.QHBoxLayout()
        layout_vsv.addLayout(tab_fname_vsv_box)
        layout_vsv.addLayout(keys_vsv)
        # ~~~~~~~~~~~ VSV END ~~~~~~~~~~~~ #


        # ~~~~~~~~~~~ DELTA START ~~~~~~~~~~~ #

        tabs_delta  =  QtWidgets.QTabWidget()
        tab1_delta  =  QtWidgets.QWidget()
        tab2_delta  =  QtWidgets.QWidget()

        frame_delta_raw  =  pg.ImageView(self, name="FrameDeltaRaw")
        frame_delta_raw.ui.roiBtn.hide()
        frame_delta_raw.ui.menuBtn.hide()
        frame_delta_raw.view.setXLink("FrameDeltaSegm")
        frame_delta_raw.view.setYLink("FrameDeltaSegm")

        frame_delta_segm  =  pg.ImageView(self, name="FrameDeltaSegm")
        frame_delta_segm.ui.roiBtn.hide()
        frame_delta_segm.ui.menuBtn.hide()
        frame_delta_segm.view.setXLink("FrameDeltaRaw")
        frame_delta_segm.view.setYLink("FrameDeltaRaw")

        frame_delta_raw_box  =  QtWidgets.QHBoxLayout()
        frame_delta_raw_box.addWidget(frame_delta_raw)

        frame_delta_segm_box  =  QtWidgets.QHBoxLayout()
        frame_delta_segm_box.addWidget(frame_delta_segm)

        tab1_delta.setLayout(frame_delta_raw_box)
        tab2_delta.setLayout(frame_delta_segm_box)

        tabs_delta.addTab(tab1_delta, "Raw Data")
        tabs_delta.addTab(tab2_delta, "Segmented")

        fname_delta_lbl  =  QtWidgets.QLabel("DELTA ", self)

        segment_delta_btn  =  QtWidgets.QPushButton("Segment", self)
        segment_delta_btn.clicked.connect(self.segment_delta)
        segment_delta_btn.setToolTip("Pop up tool to remove badly reconstructed pseudo cells")
        segment_delta_btn.setFixedSize(int(ksf_h * 130), int(ksf_w * 25))
        segment_delta_btn.setEnabled(True)

        keys_delta  =  QtWidgets.QVBoxLayout()
        keys_delta.addWidget(segment_delta_btn)
        keys_delta.addStretch()

        tab_fname_delta_box  =  QVBoxLayout()
        tab_fname_delta_box.addWidget(fname_delta_lbl)
        tab_fname_delta_box.addWidget(tabs_delta)

        layout_delta  =  QtWidgets.QHBoxLayout()
        layout_delta.addLayout(tab_fname_delta_box)
        layout_delta.addLayout(keys_delta)
        # ~~~~~~~~~~~ DELTA END ~~~~~~~~~~~ #


        # ~~~~~~~~~~~ BOTTOM WIDGETS ~~~~~~~~~~~ #
        busy_lbl  =  QtWidgets.QLabel("Ready")
        busy_lbl.setStyleSheet("color: green")

        pixsize_x_lbl  =  QtWidgets.QLabel("pix size XY =;")
        pixsize_z_lbl  =  QtWidgets.QLabel("Z step =")

        bottom_labels_box  =  QtWidgets.QHBoxLayout()
        bottom_labels_box.addWidget(busy_lbl)
        bottom_labels_box.addStretch()
        bottom_labels_box.addWidget(pixsize_x_lbl)
        bottom_labels_box.addWidget(pixsize_z_lbl)
        # ~~~~~~~~ END BOTTOM WIDGETS ~~~~~~~~~~ #

        tabs_tot   =  QtWidgets.QTabWidget()
        tab_vsv    =  QtWidgets.QWidget()
        tab_delta  =  QtWidgets.QWidget()

        tabs_tot.addTab(tab_vsv, "VSV")
        tabs_tot.addTab(tab_delta, "Delta")

        tab_vsv.setLayout(layout_vsv)
        tab_delta.setLayout(layout_delta)

        layout  =  QtWidgets.QVBoxLayout(widget)
        layout.addWidget(tabs_tot)
        layout.addLayout(bottom_labels_box)

        mycmap  =  np.fromfile("mycmap.bin", "uint16").reshape((10000, 3))      # / 255.0
        self.colors4map  =  []
        for k in range(mycmap.shape[0]):
            self.colors4map.append(mycmap[k, :])
        self.colors4map     =  self.colors4map + self.colors4map + self.colors4map + self.colors4map + self.colors4map + self.colors4map
        self.colors4map[0]  =  np.array([0, 0, 0])

        self.busy_lbl            =  busy_lbl
        self.frame_delta_raw     =  frame_delta_raw
        self.frame_delta_segm    =  frame_delta_segm
        self.frame_vsv_raw       =  frame_vsv_raw
        self.frame_vsv_segm      =  frame_vsv_segm
        self.pixsize_x_lbl       =  pixsize_x_lbl
        self.pixsize_z_lbl       =  pixsize_z_lbl
        self.fname_vsv_lbl       =  fname_vsv_lbl
        self.fname_delta_lbl     =  fname_delta_lbl

        self.thickness_thr_edt  =  thickness_thr_edt
        self.up_ratio_thr_edt   =  up_ratio_thr_edt
        self.low_ratio_thr_edt  =  low_ratio_thr_edt
        self.high_area_thr_edt  =  high_area_thr_edt

        self.soft_version  =  "TrojanHorseQuantifier_v1.0"

        self.setGeometry(800, 100, 1200, 800)
        self.setWindowTitle(self.soft_version)
        self.setWindowIcon(QtGui.QIcon('Icons/DrosophilaIcon.png'))
        self.show()

    def closeEvent(self, event):
        """Close the GUI, asking confirmation."""
        quit_msg  =  "Are you sure you want to exit the program?"
        reply     =  QtWidgets.QMessageBox()
        reply.setText(quit_msg)
        reply.setStandardButtons(QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No)
        x         =  reply.exec()

        if x == QtWidgets.QMessageBox.StandardButton.Yes:
            event.accept()
        else:
            event.ignore()

    def busy_indicator(self):
        """Write a red text (BUSY) as a label on the GUI (bottom left)."""
        self.busy_lbl.setText("Busy")
        self.busy_lbl.setStyleSheet('color: red')

    def ready_indicator(self):
        """Write a green text (READY) as a label on the GUI (bottom left)."""
        self.busy_lbl.setText("Ready")
        self.busy_lbl.setStyleSheet('color: green')

    def thickness_thr_var(self, text):
        """Set the thickness threshold value."""
        try:
            self.thickness_thr_value  =  float(text)
        except ValueError:
            pass

    def up_ratio_thr_var(self, text):
        """Set the upper ratio threshold value."""
        try:
            self.up_ratio_thr_value  =  float(text)
        except ValueError:
            pass

    def low_ratio_thr_var(self, text):
        """Set the lower ratio threshold value."""
        try:
            self.low_ratio_thr_value  =  float(text)
        except ValueError:
            pass

    def high_area_thr_var(self, text):
        """Set the uo bound area threshold value."""
        try:
            self.high_area_thr_value  =  float(text)
        except ValueError:
            pass

    def load_raw_data(self):
        """Load raw data files."""
        reload(LoadRawData)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()

        try:
            self.frame_vsv_segm.clear()
            # self.frame_vsv_flt_segm.clear()
            self.frame_delta_segm.clear()
            self.raw_vsv_fname    =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select raw data VSV file", filter='*.tif *.tiff')[0])
            self.raw_delta_fname  =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select DELTA file to couple with " + self.raw_vsv_fname[self.raw_vsv_fname.rfind('/') + 1:], filter='*.tif *.tiff')[0])
            self.fname_vsv_lbl.setText("VSV:   " + self.raw_vsv_fname[self.raw_vsv_fname.rfind('/') + 1:])
            self.fname_delta_lbl.setText("DELTA: " + self.raw_delta_fname[self.raw_delta_fname.rfind('/') + 1:])
            self.raw_data         =  LoadRawData.LoadRawData(self.raw_vsv_fname, self.raw_delta_fname)

            self.frame_vsv_raw.setImage(self.raw_data.raw_vsv)
            self.frame_delta_raw.setImage(self.raw_data.raw_delta)

            self.pixsize_x_lbl.setText("pix size XY = " + str(self.raw_data.pix_sizeX))
            self.pixsize_z_lbl.setText("pix size Z = " + str(self.raw_data.pix_sizeZ))

        except Exception:
            traceback.print_exc()

        self.ready_indicator()

    def segment_vsv(self):
        """Segment the VSV virus."""
        reload(VsvDetector)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()

        try:
            self.segm_vsv_noflt  =  VsvDetector.VsvDetector(self.raw_data.raw_vsv)
            self.frame_vsv_segm.setImage(self.segm_vsv_noflt.vsv_mask)
            self.vsv_cmap        =  pg.ColorMap(np.linspace(0, 1, self.segm_vsv_noflt.vsv_mask.max()), color=self.colors4map)
            self.frame_vsv_segm.setColorMap(self.vsv_cmap)
        except Exception:
            traceback.print_exc()

        self.ready_indicator()

    def filter_vsv(self):
        """Filter delta virus."""
        reload(VsvDetector)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()

        try:
            self.segm_vsv  =  VsvDetector.VsvFilter(self.segm_vsv_noflt.vsv_mask, self.thickness_thr_value, self.up_ratio_thr_value, self.low_ratio_thr_value, self.high_area_thr_value)
            self.frame_vsv_segm.setImage(self.segm_vsv.vsv_mask + self.segm_vsv.vsv_diff2show)
            self.vsv_cmap  =  pg.ColorMap(np.linspace(0, 1, self.segm_vsv.vsv_mask.max()), color=self.colors4map)
            self.frame_vsv_segm.setColorMap(self.vsv_cmap)
        except Exception:
            traceback.print_exc()

        self.ready_indicator()

    def segment_delta(self):
        """Segment delta virus."""
        reload(DeltaDetector)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()

        try:
            self.segm_delta  =  DeltaDetector.DeltaDetector(self.raw_data.raw_delta)
            self.frame_delta_segm.setImage(self.segm_delta.delta_spts)
            self.delta_cmap  =  pg.ColorMap(np.linspace(0, 1, self.segm_delta.delta_spts.max()), color=self.colors4map)
            self.frame_delta_segm.setColorMap(self.delta_cmap)
        except Exception:
            traceback.print_exc()

        self.ready_indicator()

    def load_analysis(self):
        """Load a previously done analysis."""
        reload(LoadRawData)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()

        try:
            self.frame_vsv_segm.clear()
            self.frame_delta_segm.clear()
            self.raw_vsv_fname    =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select raw data VSV file", filter='*.tif *.tiff')[0])
            self.raw_delta_fname  =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select DELTA file to couple with " + self.raw_vsv_fname[self.raw_vsv_fname.rfind('/') + 1:], filter='*.tif *.tiff')[0])
            self.fname_vsv_lbl.setText("VSV: " + self.raw_vsv_fname[self.raw_vsv_fname.rfind('/') + 1:])
            self.fname_delta_lbl.setText("DELTA: " + self.raw_delta_fname[self.raw_delta_fname.rfind('/') + 1:])
            self.raw_data         =  LoadRawData.LoadRawData(self.raw_vsv_fname, self.raw_delta_fname)

            self.frame_vsv_raw.setImage(self.raw_data.raw_vsv)
            self.frame_delta_raw.setImage(self.raw_data.raw_delta)

            self.pixsize_x_lbl.setText("pix size XY = " + str(self.raw_data.pix_sizeX))
            self.pixsize_z_lbl.setText("pix size Z = " + str(self.raw_data.pix_sizeZ))

            load_bff  =  np.load(self.raw_delta_fname.replace("RED", "ANALYSIS")[:-4] + "npz")
            params    =  load_bff["arr_2"]

            self.thickness_thr_edt.setText(str(params[0]))
            self.up_ratio_thr_edt.setText(str(params[1]))
            self.low_ratio_thr_edt.setText(str(params[2]))
            self.high_area_thr_edt.setText(str(params[3]))

            self.segment_vsv()
            self.filter_vsv()
            self.segment_delta()

        except Exception:
            traceback.print_exc()


    def save_analysis(self):
        """Save analysis parameters."""
        AnalysisSaver.AnalysisSaver(self.raw_delta_fname, self.segm_vsv.vsv_mask, self.segm_vsv.vsv_diff, self.segm_delta.delta_spts, self.thickness_thr_value, self.up_ratio_thr_value, self.low_ratio_thr_value, self.high_area_thr_value)

    def post_statistics(self):
        """Prepare a recap xlsx file with colocalization results."""
        analyses_folder  =  str(QtWidgets.QFileDialog.getExistingDirectory(None, "Select the Directory with the Analyses"))
        AnalysisSaver.PostStatistics(analyses_folder, self.soft_version)


def main():
    app         =  QtWidgets.QApplication(sys.argv)
    splash_pix  =  QtGui.QPixmap('Icons/TrojanHorse.jpeg')
    splash      =  QtWidgets.QSplashScreen(splash_pix, QtCore.Qt.WindowType.WindowStaysOnTopHint | QtCore.Qt.WindowType.FramelessWindowHint)
    splash.show()
    ex = MainWindow()
    QtCore.QTimer.singleShot(500, lambda: splash.finish(ex))
    sys.exit(app.exec())


def except_hook(cls, exception, traceback):
    sys.__excepthook__(cls, exception, traceback)


if __name__ == '__main__':

    main()
    sys.excepthook  =  except_hook
