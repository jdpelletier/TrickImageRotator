import os, time, sys, threading
from os import listdir
from os.path import abspath, isfile, join
from pathlib import Path
import math
import subprocess
import datetime

import numpy as np
from astropy.io import fits
import PIL.Image as PILimage

from ginga.misc import log
from ginga.qtw.QtHelp import QtGui, QtCore
from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
from ginga.util.loader import load_data

import ktl

class FitsViewer(QtGui.QMainWindow):

    def __init__(self, logger):
        super(FitsViewer, self).__init__()
        self.logger = logger

        self.cachedFiles = None
        #KTL stuff
        #Cache KTL keywords
        self.curAngle = ktl.cache('dcs', 'rotdest')
        self.trickxpos = ktl.cache('ao', 'TRKRO1XP')
        self.trickypos = ktl.cache('ao', 'TRKRO1YP')
        self.trickxsize = ktl.cache('ao', 'TRKRO1XS')
        self.trickysize = ktl.cache('ao', 'TRKRO1YS')


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
        bd = fi.get_bindings()
        bd.enable_all(True)

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

        self.readout = QtGui.QLabel("test")
        wstartscan = QtGui.QPushButton("Start Scan")
        wstartscan.clicked.connect(self.start_scan)
        wopen = QtGui.QPushButton("Open File")
        wopen.clicked.connect(self.open_file)
        wquit = QtGui.QPushButton("Quit")
        wquit.clicked.connect(self.quit)
        fi.set_callback('cursor-changed', self.motion_cb)
        hbox.addStretch(1)
        for w in (self.readout, wstartscan, wopen, wquit):
            hbox.addWidget(w, stretch=0)

        hw = QtGui.QWidget()
        hw.setLayout(hbox)
        vbox.addWidget(hw, stretch=0)

        vw = QtGui.QWidget()
        self.setCentralWidget(vw)
        vw.setLayout(vbox)
        self.dc = self.add_canvas()


    def add_canvas(self, tag=None):
        # add a canvas to the view
        my_canvas = self.fitsimage.get_canvas()
        DrawingCanvas = my_canvas.get_draw_class('rectangle')
        return DrawingCanvas

    def start_scan(self):
        hdu = fits.PrimaryHDU()
        try:
            hdu.writeto('procImage.fits')
        except OSError:
            os.remove('procImage.fits')
            hdu.writeto('procImage.fits')
        self.cachedFiles = self.walkDirectory()
        self.thread = None
        self.runThread = True
        print("scan started...")
        self.scan(1, self.cachedFiles)

    def load_file(self, filepath):
        filepath = self.processData(filepath)
        image = load_data(filepath, logger=self.logger)
        self.fitsimage.set_image(image)
#        self.setWindowTitle(filepath)
        left, right, up, down = self.getROI()
        self.box = self.dc(left, down, right, up, color='green')

        self.fitsimage.get_canvas().add(self.box, redraw=True)
        self.fitsimage.rotate(self.rotAngle())

    def open_file(self):
        res = QtGui.QFileDialog.getOpenFileName(self, "Open FITS file",
                                                ".")
        print(res)
        if isinstance(res, tuple):
            fileName = res[0]
        else:
            fileName = str(res)
        if len(fileName) != 0:
            self.load_file(fileName)

    def motion_cb(self, viewer, button, data_x, data_y):

        # Get the value under the data coordinates
        try:
            # We report the value across the pixel, even though the coords
            # change halfway across the pixel
            value = viewer.get_data(int(data_x + 0.5), int(data_y + 0.5))

        except Exception:
            value = None

        fits_x, fits_y = data_x, data_y

        # TODO ADD WCS
        # Calculate WCS RA
        # try:
        #     # NOTE: image function operates on DATA space coords
        #     image = viewer.get_image()
        #     if image is None:
        #         # No image loaded
        #         return
        #     ra_txt, dec_txt = image.pixtoradec(fits_x, fits_y,
        #                                        format='str', coords='fits')
        # except Exception as e:
        #     self.logger.warning("Bad coordinate conversion: %s" % (
        #         str(e)))
        #     ra_txt = 'BAD WCS'
        #     dec_txt = 'BAD WCS'
        #
        # text = "RA: %s  DEC: %s  X: %.2f  Y: %.2f  Value: %s" % (
        #     ra_txt, dec_txt, fits_x, fits_y, value)
        text = "X: %.2f  Y: %.2f  Value: %s" % (
             fits_x, fits_y, value)
        self.readout.setText(text)

    def quit(self, *args):
        self.logger.info("Attempting to shut down the application...")
        self.runThread = False
        try:
            self.thread.join()
        except AttributeError:
            print("Scanning never started")
        self.deleteLater()

    ##Start of image find and processing code




    #TODO:
    #add osiris FOV


    def rotAngle(self):
        cur = self.curAngle.read()
        final = 179.5 - float(cur)#TODO figure out if this angle is right
        return final

    def getROI(self):
        left = int(self.trickxpos.read())
        right = int(self.trickxpos.read()) + int(self.trickxsize.read())*5
        up = int(self.trickypos.read())
        down = int(self.trickypos.read()) + int(self.trickysize.read())*5
        print("ROI box: %d %d %d %d" %(left, right, up, down))
        return left, right, up, down

    def nightpath(self):
        # nightly = Path('/net/k1aoserver/k1aodata/nightly')
        # date = datetime.datetime.utcnow()
        # year, month, day = str(date.strftime("%y")), str(date.strftime("%m")), str(date.strftime("%d"))
        # nightly = nightly / year / month / day / 'Trick'
        return '.'

    def walkDirectory(self):
        directory = self.nightpath()
        return [abspath(join(directory, f)) for f in listdir(directory) if isfile(join(directory, f))]


    def updateFileCache(self, cachedFiles):
        updatedFileList = self.walkDirectory()
        filtered = [i for i in updatedFileList if not i in cachedFiles]
        cachedFiles = updatedFileList
        return len(filtered) > 0, filtered, cachedFiles


    def scan(self, timeout, cachedFiles):
        def __target():
            while self.runThread:
                hasNewFiles, files, self.cachedFiles = self.updateFileCache(self.cachedFiles)
                if hasNewFiles:
                    print("New Image Detected!")
                    filen = files[0]
                    self.waitForFileToBeUnlocked(filen, 1)
                    self.load_file(filen)
                time.sleep(timeout)
        self.thread = threading.Thread(target=__target)
        self.thread.daemon = True
        self.thread.start()
        return cachedFiles

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

    def processData(self, filename):
        fitsData = fits.getdata(filename, ext=0)
        header = fits.getheader(filename)
        mask = fits.getdata('BadPix_1014Hz.fits', ext=0)
        maskedData = np.multiply(fitsData, mask)
        background = np.median(maskedData)
        return self.writeFits(header, maskedData - background)

    def writeFits(self, headerinfo, image_data):
        hdu = fits.PrimaryHDU()
        hdu.data = image_data
        hdu.header = headerinfo
        filename = 'procImage.fits'
        try:
            hdul = fits.HDUList([hdu])
            hdul.writeto(filename)
            hdul.close()
        except OSError:
            os.remove(filename)
            hdul = fits.HDUList([hdu])
            hdul.writeto(filename)
            hdul.close()
        return filename

def main():
    ##Write dummy file so walkDirectory caches it in the beginning

    app = QtGui.QApplication([])

    # ginga needs a logger.
    # If you don't want to log anything you can create a null logger by
    # using null=True in this call instead of log_stderr=True
    logger = log.get_logger("example1", log_stderr=True, level=40)

    w = FitsViewer(logger)
    w.resize(700, 800)
    w.show()
    app.setActiveWindow(w)
    w.raise_()
    w.activateWindow()
    app.exec_()

if __name__ == "__main__":
    main()
