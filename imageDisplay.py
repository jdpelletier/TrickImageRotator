import os, time, sys
from os import listdir
from os.path import abspath, isfile, join
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
from matplotlib import colors
import PIL.Image as PILimage
import math

#import ktl

plt.style.use(astropy_mpl_style)

#curAngle = ktl.cache('dcs', 'ROTPOSN')
#parAngle = ktl.cache('dcs', 'PARANG')

#def rotAngle(current, para):
#    cur = current.read()
#    par = para.read()
#    return par - cur

def walkDirectory():
    directory = '.'
    return [abspath(join(directory, f)) for f in listdir(directory) if isfile(join(directory, f))]


def updateFileCache(cachedFiles):
    updatedFileList = walkDirectory()
    filtered = [i for i in updatedFileList if not i in cachedFiles]
    cachedFiles = updatedFileList
    return len(filtered) > 0, filtered, cachedFiles


def scan(timeout, cachedFiles):
    hasNewFiles, files, cachedFiles = updateFileCache(cachedFiles)
    if hasNewFiles:
        print("New File Detected!")
        filen = files[0]
        waitForFileToBeUnlocked(filen, 1)
        fitsData = fits.getdata(filen, ext=0)
        header = fits.getheader(filen)
        newdata = rotate_clip(fitsData, 45, 942, 747)
        writeFits(header, newdata)
        print("File closed")
    time.sleep(timeout)
    return cachedFiles

def fileIsCurrentlyLocked(filepath):
    locked = None
    hdulist = None
    file_object = None
    if os.path.exists(filepath):
        try:
            print("Trying to open %s." % filepath)

            hdulist = fits.open(filepath)

            file_object = np.sum([1 for hdu in hdulist if type(hdu) in
                    	[fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU]
                    	and hdu.data is not None])
            if file_object:
                print("%s is not locked." % filepath)
                locked = False

        except TypeError:
            print("File is locked (unable to open in append mode)")
            locked = True

        finally:
            if file_object:
                hdulist.close()
                print("%s closed." % filepath)

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
    '''
    #numpy rot
    yi, xi = np.mgrid[0:ht, 0:wd]
    xi -= rotctr_x
    yi -= rotctr_y
    cos_t = np.cos(np.radians(theta_deg))
    sin_t = np.sin(np.radians(theta_deg))

    if have_numexpr:
        ap = ne.evaluate("(xi * cos_t) - (yi * sin_t) + rotctr_x")
        bp = ne.evaluate("(xi * sin_t) + (yi * cos_t) + rotctr_y")
    else:
        ap = (xi * cos_t) - (yi * sin_t) + rotctr_x
        bp = (xi * sin_t) + (yi * cos_t) + rotctr_y

    #ap = np.rint(ap).clip(0, wd-1).astype(np.int)
    #bp = np.rint(bp).clip(0, ht-1).astype(np.int)
    # Optomizations to reuse existing intermediate arrays
    np.rint(ap, out=ap)
    ap = ap.astype(np.int, copy=False)
    ap.clip(0, wd - 1, out=ap)
    np.rint(bp, out=bp)
    bp = bp.astype(np.int, copy=False)
    bp.clip(0, ht - 1, out=bp)

    if out is not None:
        out[:, :, ...] = data_np[bp, ap]
        newdata = out
    else:
        newdata = data_np[bp, ap]
        new_ht, new_wd = newdata.shape[:2]

        assert (wd == new_wd) and (ht == new_ht), \
            Exception("rotated cutout is %dx%d original=%dx%d" % (
                new_wd, new_ht, wd, ht))
    #end np rot
    '''
    return newdata

def writeFits(headerinfo, image_data):
    hdu = fits.PrimaryHDU()
    hdu.data = image_data
    hdu.header = headerinfo
    try:
        hdu.writeto('rotatedImage.fits')
    except OSError:
        os.remove('rotatedImage.fits')
        hdu.writeto('rotatedImage.fits')


#####

def main():
    cachedFiles = None
    cachedFiles = walkDirectory();
    print("Scan started...")
    try:
        while True:
            cachedFiles = scan(1, cachedFiles)
    except KeyboardInterrupt:
            stop_scan()

if __name__ == "__main__":
    main()
