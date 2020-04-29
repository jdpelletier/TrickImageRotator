import os, time, sys
from os import listdir
from os.path import abspath, isfile, join
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
from matplotlib import colors

plt.style.use(astropy_mpl_style)



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
        waitForFileToBeUnlocked(filen, 1);
        plotIm(filename=filen)
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


def plotIm(filename):
    image_data = fits.getdata(filename, ext=0)
    plt.figure()
    plt.imshow(image_data, cmap='gray', norm=colors.LogNorm())
    plt.colorbar()
    plt.show()

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
