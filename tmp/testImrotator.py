import os
import argparse
from astropy.io import fits
import math
import numpy as np
import PIL.Image as PILimage
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(description="Rotates image",
                         usage="imageRotator.py file angle")

parser.add_argument("file", help="file to load")
parser.add_argument("angle", type=int, help="angle to rotate")

args = parser.parse_args()


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

def plotIm(image_data):
    plt.figure()
    plt.imshow(image_data, cmap='gray')
    plt.colorbar()
    plt.show()

def writeFits(headerinfo, image_data):
    hdu = fits.PrimaryHDU()
    hdu.data = image_data
    hdu.header = headerinfo
    try:
        hdu.writeto('rotatedImage.fits')
    except OSError:
        os.remove('rotatedImage.fits')
        hdu.writeto('rotatedImage.fits')

def main():
    f = fits.open(args.file)
    fitsData = fits.getdata(args.file, ext=0)
    header = fits.getheader(args.file)
    newdata = rotate_clip(fitsData, args.angle, 942, 747)
    writeFits(header, fitsData)
    plotIm(fitsData)
    plotIm(newdata)


if __name__=='__main__':
    main()
