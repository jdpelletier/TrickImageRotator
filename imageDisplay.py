import os, time, sys
from os import listdir
from os.path import abspath, isfile, join
from pathlib import Path
import math
import subprocess
import datetime

import numpy as np
from astropy.io import fits
import PIL.Image as PILimage


import ktl

#Cache KTL keywords
curAngle = ktl.cache('dcs', 'rotdest')
trickxpos = ktl.cache('ao', 'TRKRO1XP')
trickypos = ktl.cache('ao', 'TRKRO1YP')
trickxsize = ktl.cache('ao', 'TRKRO1XS')
trickysize = ktl.cache('ao', 'TRKRO1YS')


#TODO:
#add osiris FOV


def rotAngle():
    cur = curAngle.read()
    final = 179.5 - float(cur)#TODO figure out if this angle is right
    return final

def getROI():
    left = int(trickxpos.read())
    right = int(trickxpos.read()) + int(trickxsize.read())*5
    up = int(trickypos.read())
    down = int(trickypos.read()) + int(trickysize.read())*5
    return left, right, up, down

def nightpath():
    nightly = Path('/net/k1aoserver/k1aodata/nightly')
    date = datetime.datetime.utcnow()
    year, month, day = str(date.strftime("%y")), str(date.strftime("%m")), str(date.strftime("%d"))
    nightly = nightly / year / month / day / 'Trick'
    return nightly

def walkDirectory():
    directory = nightpath()
    return [abspath(join(directory, f)) for f in listdir(directory) if isfile(join(directory, f))]


def updateFileCache(cachedFiles):
    updatedFileList = walkDirectory()
    filtered = [i for i in updatedFileList if not i in cachedFiles]
    cachedFiles = updatedFileList
    return len(filtered) > 0, filtered, cachedFiles


def scan(timeout, cachedFiles):
    hasNewFiles, files, cachedFiles = updateFileCache(cachedFiles)
    if hasNewFiles:
        print("New Image Detected!")
        filen = files[0]
        waitForFileToBeUnlocked(filen, 1)
        fitsData = fits.getdata(filen, ext=0)
        header = fits.getheader(filen)
        procdata = processImage(fitsData)
        procdata, max = buildROIBox(procdata)
        print("Rotating image %f degrees" % rotAngle())
        procdata = rotate_clip(procdata, rotAngle(), 942, 747)
        displayFits(header, procdata, max)
    time.sleep(timeout)
    return cachedFiles

def fileIsCurrentlyLocked(filepath):
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

def waitForFileToBeUnlocked(filename, wait_time):
    # if the file doesn't exist, wait wait_time seconds and try again until it's found
    while not os.path.exists(filename):
        print("%s hasn't arrived. Waiting %s seconds." % (filename, wait_time))
        time.sleep(wait_time)

    # if the file exists but locked, wait wait_time seconds and check
    # again until it's no longer locked by another process
    while fileIsCurrentlyLocked(filename):
        print("%s is currently in use. Waiting %s seconds." % (filename, wait_time))
        time.sleep(wait_time)



def stop_scan():
    print("Shutting down...")
    sys.exit()


def processImage(fitsdata):
    mask = fits.getdata('BadPix_1014Hz.fits', ext=0)
    maskedData = np.multiply(fitsdata, mask)
    background = np.median(maskedData)
    return maskedData - background

def rotate_clip(data_np, theta_deg, rotctr_x=None, rotctr_y=None,
                out=None, logger=None):
    """
    Rotate numpy array `data_np` by `theta_deg` around rotation center
    (rotctr_x, rotctr_y).  If the rotation center is omitted it defaults
    to the center of the array.
    No adjustment is done to the data array beforehand, so the result will
    be clipped according to the size of the array (the output array will be
    the same size as the input array).
    """

    # If there is no rotation, then we are done
    if math.fmod(theta_deg, 360.0) == 0.0:
        return data_np

    ht, wd = data_np.shape[:2]
    dtype = data_np.dtype

    if rotctr_x is None:
        rotctr_x = wd // 2
    if rotctr_y is None:
        rotctr_y = ht // 2

    #pil rot
    img = PILimage.fromarray(data_np)
    img_rot = img.rotate(theta_deg, resample=False, expand=False,
                         center=(rotctr_x, rotctr_y))
    newdata = np.array(img_rot, dtype=data_np.dtype)
    new_ht, new_wd = newdata.shape[:2]
    assert (wd == new_wd) and (ht == new_ht), \
        Exception("rotated cutout is %dx%d original=%dx%d" % (
            new_wd, new_ht, wd, ht))
    return newdata

def buildROIBox(data):
    left, right, up, down = getROI()
    max = np.amax(data) + 1.0
    data[left][up] = max
    data[right][up] = max
    data[right][down] = max
    data[left][down] = max
    return data, max

def findBox(data, max):
    result = np.where(data == max)
    print(result)
    ty = result[0][0]
    tx = result[1][0]
    ly = result[0][1]
    lx = result[1][1]
    ry = result[0][2]
    rx = result[1][2]
    by = result[0][3]
    bx = result[1][3]
    return tx, ty, lx, ly, rx, ry, bx, by

def writeFits(headerinfo, image_data):
    hdu = fits.PrimaryHDU()
    hdu.data = image_data
    hdu.header = headerinfo
    filename = 'rotatedImage.fits'
    try:
        hdu.writeto(filename)
    except OSError:
        os.remove(filename)
        hdu.writeto(filename)
    return filename

def displayFits(headerinfo, newdata, max):
    pgrep = subprocess.Popen("pgrep ds9_80", stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
    (output, err) = pgrep.communicate()
    if output != '':
        subprocess.Popen("pkill ds9_80", stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
    tx, ty, lx, ly, rx, ry, bx, by = findBox(newdata, max)
    filename = writeFits(headerinfo, newdata)
    command = "ds9_80 %s -scale SQRT -zoom TO FIT -regions command 'line 240 1970 240 1790' " \
    "-regions command 'line 60 1790 240 1790' -regions command 'line 210 1940 240 1970' " \
    "-regions command 'line 240 1970 270 1940' -regions command 'line 90 1820 60 1790' " \
    "-regions command 'line 90 1760 60 1790' -regions command 'text 290 1880 #text=\"N\"' " \
    "-regions command 'line %d %d %d %d #color = red' -regions command 'line %d %d %d %d #color = red' " \
    "-regions command 'line %d %d %d %d #color = red' -regions command 'line %d %d %d %d #color = red' " \
    "-regions command 'text 150 1710 #text=\"E\"' -regions command 'line 924 100 1124 100' " \
    "-regions command 'text 1024 50 #text= \"10 as\" font=\"bold\"'" % (filename, lx, ly, tx, ty, tx, ty, rx, ry, rx, ry, bx, by, bx, by, lx, ly)
    subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)

def main():
    ##Write dummy file so walkDirectory caches it in the beginning
    hdu = fits.PrimaryHDU()
    try:
        hdu.writeto('rotatedImage.fits')
    except OSError:
        os.remove('rotatedImage.fits')
        hdu.writeto('rotatedImage.fits')
    cachedFiles = None
    cachedFiles = walkDirectory()
    print("Scan started...")
    try:
        while True:
            cachedFiles = scan(1, cachedFiles)
    except KeyboardInterrupt:
            stop_scan()

if __name__ == "__main__":
    main()
