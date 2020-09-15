import os, time, sys, threading
from os import listdir
from os.path import abspath, isfile, join
from pathlib import Path
import math
import subprocess
import datetime

import numpy as np
from astropy.io import fits
from astropy import wcs
import PIL.Image as PILimage

from ginga import Bindings
from ginga.misc import log
from ginga.qtw.QtHelp import QtGui, QtCore
from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
from ginga.util.loader import load_data

import ktl

class ScannerSignals(QtCore.QObject):
    load = QtCore.Signal(object)

class Scanner(QtCore.QRunnable):
    '''
    Scanner thread
    '''
    def __init__(self, fn, *args, **kwargs):
        super(Scanner, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = ScannerSignals()
        self.kwargs['file_callback'] = self.signals.load

        # Add the callback to our kwargs
    @QtCore.Slot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        self.fn(*self.args, **self.kwargs)


class FitsViewer(QtGui.QMainWindow):

    def __init__(self, logger):
        super(FitsViewer, self).__init__()
        self.logger = logger

        self.cachedFiles = None
        self.scanning = False
        #KTL stuff
        #Cache KTL keywords
        self.trickxpos = ktl.cache('ao', 'TRKRO1XP')
        self.trickypos = ktl.cache('ao', 'TRKRO1YP')
        self.trickxsize = ktl.cache('ao', 'TRKRO1XS')
        self.trickysize = ktl.cache('ao', 'TRKRO1YS')

        self.threadpool = QtCore.QThreadPool()

        # create the ginga viewer and configure it
        fi = CanvasView(self.logger, render='widget')
        fi.enable_autocuts('on')
        fi.set_autocut_params('zscale')
        fi.enable_autozoom('on')
        # fi.set_callback('drag-drop', self.drop_file)
        fi.set_bg(0.2, 0.2, 0.2)
        fi.ui_set_active(True)
        self.fitsimage = fi

        # enable some user interaction
        self.bd = fi.get_bindings()
        self.bd.enable_all(True)

        w = fi.get_widget()
        w.resize(512, 512)

        # add scrollbar interface around this viewer
        si = ScrolledView(fi)

        vbox = QtGui.QVBoxLayout()
        vbox.setContentsMargins(QtCore.QMargins(2, 2, 2, 2))
        vbox.setSpacing(1)
        vbox.addWidget(si, stretch=1)

        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(QtCore.QMargins(4, 2, 4, 2))

        self.readout = QtGui.QLabel("")

        hbox.addStretch(1)
        hbox.addWidget(self.readout, stretch = 0)

        hw = QtGui.QWidget()
        hw.setLayout(hbox)
        vbox.addWidget(hw, stretch=0)

        hbox2 = QtGui.QHBoxLayout()
        hbox2.setContentsMargins(QtCore.QMargins(4, 2, 4, 2))
        self.wstartscan = QtGui.QPushButton("Start Scan")
        self.wstartscan.clicked.connect(self.start_scan)
        self.wstartscan.setEnabled(True)
        self.wstopscan = QtGui.QPushButton("Stop Scan")
        self.wstopscan.clicked.connect(self.stop_scan)
        self.wstopscan.setEnabled(False)
        self.wcut = QtGui.QComboBox()
        for name in fi.get_autocut_methods():
            self.wcut.addItem(name)
        self.wcut.currentIndexChanged.connect(self.cut_change)
        self.wcolor = QtGui.QComboBox()
        for name in fi.get_color_algorithms():
            self.wcolor.addItem(name)
        self.wcolor.currentIndexChanged.connect(self.color_change)
        wopen = QtGui.QPushButton("Open File")
        wopen.clicked.connect(self.open_file)
        wquit = QtGui.QPushButton("Quit")
        wquit.clicked.connect(self.quit)
        fi.set_callback('cursor-changed', self.motion_cb)
        fi.add_callback('cursor-down', self.btndown)
        hbox2.addStretch(1)
        for w in (self.wstartscan, self.wstopscan, self.wcut, self.wcolor, wopen, wquit):
            hbox2.addWidget(w, stretch=0)

        hw2 = QtGui.QWidget()
        hw2.setLayout(hbox2)
        vbox.addWidget(hw2, stretch=0)

        vw = QtGui.QWidget()
        self.setCentralWidget(vw)
        vw.setLayout(vbox)
        self.recdc, self.compdc = self.add_canvas()
        self.boxtag = "roi-box"


    def add_canvas(self, tag=None):
        # add a canvas to the view
        my_canvas = self.fitsimage.get_canvas()
        RecCanvas = my_canvas.get_draw_class('rectangle')
        CompCanvas = my_canvas.get_draw_class('compass')
        return RecCanvas, CompCanvas

    def start_scan(self):
        self.wstartscan.setEnabled(False)
        self.wstopscan.setEnabled(True)
        self.scanning = True
        hdu = fits.PrimaryHDU()
        try:
            hdu.writeto('procImage.fits')
        except OSError:
            os.remove('procImage.fits')
            hdu.writeto('procImage.fits')
        self.cachedFiles = self.walkDirectory()
        print("scan started...")
        scanner = Scanner(self.scan)
        scanner.signals.load.connect(self.load_file)
        self.threadpool.start(scanner)

    def stop_scan(self):
        self.wstartscan.setEnabled(True)
        self.wstopscan.setEnabled(False)
        self.scanning = False
        print('Scanning stopped.')

    def load_file(self, filepath):
        filepath = self.processData(filepath)
        image = load_data(filepath, logger=self.logger)
        self.fitsimage.set_image(image)
        # self.setWindowTitle(filepath)
        left, right, up, down = self.getROI()
        try:
            self.fitsimage.get_canvas().get_object_by_tag(self.boxtag)
            self.fitsimage.get_canvas().delete_object_by_tag(self.boxtag)
            self.box = self.recdc(left, down, right, up, color='green')
            self.fitsimage.get_canvas().add(self.box, tag=self.boxtag, redraw=True)
        except KeyError:
            self.box = self.recdc(left, down, right, up, color='green')
            self.fitsimage.get_canvas().add(self.box, tag=self.boxtag, redraw=True)
        width, height = image.get_size()
        data_x, data_y = width / 2.0, height / 2.0
        # x, y = self.fitsimage.get_canvas_xy(data_x, data_y)
        radius = float(max(width, height)) / 20
        self.fitsimage.get_canvas().add(self.compdc(data_x, data_y, radius, color='skyblue',
                                       fontsize=8))
        self.bd._orient(self.fitsimage, righthand=False, msg=True)

    def open_file(self):
        res = QtGui.QFileDialog.getOpenFileName(self, "Open FITS file",
                                                str(self.nightpath()))
        print(res)
        if isinstance(res, tuple):
            fileName = res[0]
        else:
            fileName = str(res)
        if len(fileName) != 0:
            self.load_file(fileName)

    def cut_change(self):
        self.fitsimage.set_autocut_params(self.wcut.currentText())

    def color_change(self):
        self.fitsimage.set_color_algorithm(self.wcolor.currentText())

    def motion_cb(self, viewer, button, data_x, data_y):

        # Get the value under the data coordinates
        try:
            # We report the value across the pixel, even though the coords
            # change halfway across the pixel
            value = viewer.get_data(int(data_x + 0.5), int(data_y + 0.5))

        except Exception:
            value = None

        fits_x, fits_y = data_x, data_y

        # Calculate WCS RA
        try:
            # NOTE: image function operates on DATA space coords
            image = viewer.get_image()
            if image is None:
                # No image loaded
                return
            ra_txt, dec_txt = image.pixtoradec(fits_x, fits_y,
                                               format='str', coords='fits')
        except Exception as e:
            self.logger.warning("Bad coordinate conversion: %s" % (
                str(e)))
            ra_txt = 'BAD WCS'
            dec_txt = 'BAD WCS'

        text = "RA: %s  DEC: %s  X: %.2f  Y: %.2f  Value: %s" % (
            ra_txt, dec_txt, fits_x, fits_y, value)
        self.readout.setText(text)

    def btndown(self, canvas, event, data_x, data_y):
        self.fitsimage.set_pan(data_x, data_y)

    def quit(self, *args):
        self.logger.info("Attempting to shut down the application...")
        self.threadpool = False
        self.deleteLater()

    ##Start of image find and processing code

    def scan(self, file_callback):
        while self.scanning:
            hasNewFiles, files, self.cachedFiles = self.updateFileCache(self.cachedFiles)
            if hasNewFiles and ('.fits' in files[0] or '.fits.gz' in files[0]):
                print("New Image Detected!")
                filen = files[0]
                self.waitForFileToBeUnlocked(filen, 1)
                file_callback.emit(filen)
            time.sleep(1)

    def walkDirectory(self):
        directory = self.nightpath()
        return [abspath(join(directory, f)) for f in listdir(directory) if isfile(join(directory, f))]

    def updateFileCache(self, cachedFiles):
        updatedFileList = self.walkDirectory()
        filtered = [i for i in updatedFileList if not i in cachedFiles]
        cachedFiles = updatedFileList
        return len(filtered) > 0, filtered, cachedFiles

    def fileIsCurrentlyLocked(self, filepath):
        locked = None
        hdulist = None
        file_object = None
        if os.path.exists(filepath):
            try:
                print("Trying to open %s." % filepath)
                #time.sleep(15) #place holder if OSError catch doesn't work
                hdulist = fits.open(filepath)

                file_object = np.sum([1 for hdu in hdulist if type(hdu) in
                        	[fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU]
                        	and hdu.data is not None])
                if file_object:
                    locked = False

            except TypeError:
                locked = True

            except OSError:
                locked = True

            finally:
                if file_object:
                    hdulist.close()

        else:
            print("%s not found." % filepath)

        return locked


    #    Checks if the files are ready.
    #    For a file to be ready it must exist and can be opened in append mode.

    def waitForFileToBeUnlocked(self, filename, wait_time):
        # if the file doesn't exist, wait wait_time seconds and try again until it's found
        while not os.path.exists(filename):
            print("%s hasn't arrived. Waiting %s seconds." % (filename, wait_time))
            time.sleep(wait_time)

        # if the file exists but locked, wait wait_time seconds and check
        # again until it's no longer locked by another process
        while self.fileIsCurrentlyLocked(filename):
            print("%s is currently in use. Waiting %s seconds." % (filename, wait_time))
            time.sleep(wait_time)

    def getROI(self):
        left = int(self.trickxpos.read()) - int(self.trickxsize.read())*3
        right = int(self.trickxpos.read()) + int(self.trickxsize.read())*3
        up = int(self.trickypos.read()) - int(self.trickysize.read())*3
        down = int(self.trickypos.read()) + int(self.trickysize.read())*3
        print("ROI box: %d %d %d %d" %(left, right, up, down))
        return left, right, up, down

    def nightpath(self):
        nightly = Path('/net/k1aoserver/k1aodata/nightly')
        date = datetime.datetime.utcnow()
        year, month, day = str(date.strftime("%y")), str(date.strftime("%m")), str(date.strftime("%d"))
        nightly = nightly / year / month / day / 'Trick'
        return nightly

    def processData(self, filename):
        header, fitsData, filter = self.addWcs(filename)
        mask = fits.getdata('BadPix_1014Hz.fits', ext=0)
        maskedData = np.multiply(fitsData, mask)
        if filter == 'H':
            background = fits.getdata('/kroot/rel/ao/qfix/data/Trick/H_sky.fits')
        else:
            background = fits.getdata('/kroot/rel/ao/qfix/data/Trick/ks_sky.fits')
        return self.writeFits(header, maskedData - background)

    def addWcs(self, filen):
        w = wcs.WCS(naxis=2)
        fitsData = fits.getdata(filen, ext=0)
        header = fits.getheader(filen)
        ht, wd = fitsData.shape[:2]
        y = ht//2
        x = wd//2
        ra = float(header['RA'])
        dec = float(header['DEC'])
        rot = float(header['ROTPOSN'])
        filter = header['TRFWNAME']
        w.wcs.crpix = [y, x]
        w.wcs.cdelt = np.array([-0.05, 0.05])
        w.wcs.crota = np.array([0.05, rot])
        w.wcs.crval = [ra, dec]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        pixcrd = np.array([[0, 0], [24, 38], [45, 98]], dtype=np.float64)
        world = w.wcs_pix2world(pixcrd, 0)
        # Convert the same coordinates back to pixel coordinates.
        pixcrd2 = w.wcs_world2pix(world, 0)
        # These should be the same as the original pixel coordinates, modulo
        # some floating-point error.
        assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6
        # Now, write out the WCS object as a FITS header
        header = w.to_header()
        return header, fitsData, filter

    def writeFits(self, headerinfo, image_data):
        hdu = fits.PrimaryHDU(header=headerinfo, data=image_data)
        filename = 'procImage.fits'
        try:
            hdu.writeto(filename)
        except OSError:
            os.remove(filename)
            hdu.writeto(filename)
        return filename

def main():
    ##Write dummy file so walkDirectory caches it in the beginning

    app = QtGui.QApplication([])

    # ginga needs a logger.
    # If you don't want to log anything you can create a null logger by
    # using null=True in this call instead of log_stderr=True
    logger = log.get_logger("example1", log_stderr=True, level=40)

    w = FitsViewer(logger)
    w.resize(600, 700)
    w.show()
    app.setActiveWindow(w)
    w.raise_()
    w.activateWindow()
    app.exec_()

if __name__ == "__main__":
    main()
