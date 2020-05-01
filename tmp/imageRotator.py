import argparse
from astropy.io import fits
import math
import numpy as np
from astropy import wcs
import PIL.Image as PILimage
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(description="Rotates image",
                         usage="imageRotator.py file")

parser.add_argument("file", help="file to load")


args = parser.parse_args()


def _array_converter(func, sky, *args, ra_dec_order=False):
    """
    A helper function to support reading either a pair of arrays
    or a single Nx2 array.
    """

    def _return_list_of_arrays(axes, origin):
        if any([x.size == 0 for x in axes]):
            return axes

        try:
            axes = np.broadcast_arrays(*axes)
        except ValueError:
            raise ValueError(
                "Coordinate arrays are not broadcastable to each other")

        xy = np.hstack([x.reshape((x.size, 1)) for x in axes])

        if ra_dec_order and sky == 'input':
            xy = wcs.WCS._denormalize_sky(xy)
        output = func(xy, origin)
        if ra_dec_order and sky == 'output':
            output = wcs.WCS._normalize_sky(output)
            return (output[:, 0].reshape(axes[0].shape),
                    output[:, 1].reshape(axes[0].shape))
        return [output[:, i].reshape(axes[0].shape)
                for i in range(output.shape[1])]

    def _return_single_array(xy, origin):
        if 0 in xy.shape:
            return xy
        if ra_dec_order and sky == 'input':
            xy = wcs.WCS._denormalize_sky(xy)
        result = func(xy, origin)
        if ra_dec_order and sky == 'output':
            result = wcs.WCS._normalize_sky(result)
        return result

    xy, origin = args
    xy = np.asarray(xy)
    origin = int(origin)

    if xy.shape == () or len(xy.shape) == 1:
        return _return_list_of_arrays([xy], origin)
    return _return_single_array(xy, origin)



def wcsradectopix(ra_deg, dec_deg, w, coords='data', naxispath=None):

        if coords == 'data':
            origin = 0
        else:
            origin = 1

        args = [ra_deg, dec_deg]
        if naxispath:
            args += [0] * len(naxispath)
        skycrd = np.array([args], np.float_)

        pix = w.wcs_world2pix(skycrd, origin)

        x = float(pix[0, 0])
        y = float(pix[0, 1])
        return (x, y)

def wcspixtoradec(w, idxs, coords='data'):

        if coords == 'data':
            origin = 0
        else:
            origin = 1
        pixcrd = np.array([idxs], np.float_)
        # sky = wcs.wcs_pix2sky(pixcrd, origin)
        # sky = wcs.all_pix2sky(pixcrd, origin)
        # astropy only?
        #sky = wcs.WCS.all_pix2world(pixcrd, origin)

        sky = _array_converter(w._all_pix2world, 'output', pixcrd, origin)
        ra_deg = float(sky[0, 0])
        dec_deg = float(sky[0, 1])

        return ra_deg, dec_deg

def dmsToDeg(sign, deg, min, sec):
    """Convert dec sign, degrees, minutes, seconds into a signed angle in
    degrees."""
    return sign * (deg + min * degPerDmsMin + sec * degPerDmsSec)



def decTimeToDeg(sign_sym, deg, min, sec):
    """Convert dec sign, degrees, minutes, seconds into a signed angle in
    degrees.
    ``sign_sym`` may represent negative as either '-' or numeric -1."""
    if sign_sym == -1 or sign_sym == '-':
        sign = -1
    else:
        sign = 1
    return dmsToDeg(sign, deg, min, sec)

def hmsStrToDeg(ra):
    """Convert a string representation of RA into a float in degrees."""
    hour, min, sec = ra.split(':')
    ra_deg = hmsToDeg(int(hour), int(min), float(sec))
    return ra_deg


def dmsStrToDeg(dec):
    """Convert a string representation of DEC into a float in degrees."""
    sign_deg, min, sec = dec.split(':')
    sign = sign_deg[0:1]
    if sign not in ('+', '-'):
        sign = '+'
        deg = sign_deg
    else:
        deg = sign_deg[1:]
    dec_deg = decTimeToDeg(sign, int(deg), int(min), float(sec))
    return dec_deg


def lon_to_deg(lon):
    """Convert longitude to degrees."""
    if isinstance(lon, str) and (':' in lon):
        # TODO: handle other coordinate systems
        lon_deg = hmsStrToDeg(lon)
    else:
        lon_deg = float(lon)
    return lon_deg

def lat_to_deg(lat):
    """Convert latitude to degrees."""
    if isinstance(lat, str) and (':' in lat):
        # TODO: handle other coordinate systems
        lat_deg = dmsStrToDeg(lat)
    else:
        lat_deg = float(lat)
    return lat_deg


def radectopix(ra_deg, dec_deg, w, format='deg', coords='data'):
    if format != 'deg':
        # convert coordinates to degrees
        ra_deg = wcs.WCS.lon_to_deg(ra_deg)
        dec_deg = wcs.WCS.lat_to_deg(dec_deg)
    return wcsradectopix(ra_deg, dec_deg, w, coords=coords,
                               naxispath=None)


def degToDms(dec, isLatitude=True):
    """Convert the dec, in degrees, to an (sign,D,M,S) tuple.
    D and M are integer, and sign and S are float.
    """
    if isLatitude:
        assert dec <= 90, WCSError("DEC (%f) > 90.0" % (dec))
        assert dec >= -90, WCSError("DEC (%f) < -90.0" % (dec))

    if dec < 0.0:
        sign = -1.0
    else:
        sign = 1.0
    dec = dec * sign

    #mnt = (dec % 1.0) * 60.0
    #sec = (dec % (1.0/60.0)) * 3600.0
    # this calculation with return values produces conversion problem.
    # e.g. dec +311600.00 -> 31.2666666667 degree
    # deg=31 min=15 sec=60 instead deg=31 min=16 sec=0.0
    # bug fixed
    mnt, sec = divmod(dec * 3600, 60)
    deg, mnt = divmod(mnt, 60)

    return (int(sign), int(deg), int(mnt), sec)

def degToHms(ra):
    """Converts the ra (in degrees) to HMS three tuple.
    H and M are in integer and the S part is in float.
    """
    assert (ra >= 0.0), WCSError("RA (%f) is negative" % (ra))
    assert ra < 360.0, WCSError("RA (%f) > 360.0" % (ra))
    rah = ra / degPerHMSHour
    ramin = (ra % degPerHMSHour) * HMSMinPerDeg
    rasec = (ra % degPerHMSMin) * HMSSecPerDeg
    return (int(rah), int(ramin), rasec)

def deg2fmt(ra_deg, dec_deg, format):
    """Format coordinates."""

    rhr, rmn, rsec = degToHms(ra_deg)
    dsgn, ddeg, dmn, dsec = degToDms(dec_deg)

    if format == 'hms':
        return rhr, rmn, rsec, dsgn, ddeg, dmn, dsec

    elif format == 'str':
        #ra_txt = '%02d:%02d:%06.3f' % (rhr, rmn, rsec)
        ra_txt = '%d:%02d:%06.3f' % (rhr, rmn, rsec)
        if dsgn < 0:
            dsgn = '-'
        else:
            dsgn = '+'
        #dec_txt = '%s%02d:%02d:%05.2f' % (dsgn, ddeg, dmn, dsec)
        dec_txt = '%s%d:%02d:%05.2f' % (dsgn, ddeg, dmn, dsec)
        return ra_txt, dec_txt


def pixtoradec(x, y, w, format='deg', coords='data'):
    args = [x, y]
    ra_deg, dec_deg = wcspixtoradec(w=w, idxs = args, coords=coords)

    if format == 'deg':
        return ra_deg, dec_deg
    return wcs.WCS.deg2fmt(ra_deg, dec_deg, format)

def add_offset_radec(ra_deg, dec_deg, delta_deg_ra, delta_deg_dec):
    """
    Algorithm to compute RA/Dec from RA/Dec base position plus tangent
    plane offsets.
    """
    # To radians
    x = math.radians(delta_deg_ra)
    y = math.radians(delta_deg_dec)
    raz = math.radians(ra_deg)
    decz = math.radians(dec_deg)

    sdecz = math.sin(decz)
    cdecz = math.cos(decz)

    d = cdecz - y * sdecz

    ra2 = math.atan2(x, d) + raz
    # Normalize ra into the range 0 to 2*pi
    twopi = math.pi * 2
    ra2 = math.fmod(ra2, twopi)
    if ra2 < 0.0:
        ra2 += twopi
    dec2 = math.atan2(sdecz + y * cdecz, math.sqrt(x * x + d * d))

    # back to degrees
    ra2_deg = math.degrees(ra2)
    dec2_deg = math.degrees(dec2)

    return (ra2_deg, dec2_deg)

def add_offset_xy(image, x, y, delta_deg_x, delta_deg_y, w):
    # calculate ra/dec of x,y pixel
    ra_deg, dec_deg = pixtoradec(x, y, w)

    # add offsets
    ra2_deg, dec2_deg = add_offset_radec(ra_deg, dec_deg,
                                         delta_deg_x, delta_deg_y)

    # then back to new pixel coords
    x2, y2 = radectopix(ra2_deg, dec2_deg, w)

    return (x2, y2)


def calc_compass(image, x, y, len_deg_e, len_deg_n, w):

    # Get east and north coordinates
    xe, ye = add_offset_xy(image, x, y, len_deg_e, 0.0, w)
    xn, yn = add_offset_xy(image, x, y, 0.0, len_deg_n, w)
    return (x, y, xn, yn, xe, ye)


def calc_compass_radius(image, x, y, radius_px, w):
    xe, ye = add_offset_xy(image, x, y, 1.0, 0.0, w)
    xn, yn = add_offset_xy(image, x, y, 0.0, 1.0, w)

    # now calculate the length in pixels of those arcs
    # (planar geometry is good enough here)
    px_per_deg_e = math.sqrt(math.fabs(ye - y)**2 + math.fabs(xe - x)**2)
    px_per_deg_n = math.sqrt(math.fabs(yn - y)**2 + math.fabs(xn - x)**2)

    # now calculate the arm length in degrees for each arm
    # (this produces same-length arms)
    len_deg_e = radius_px / px_per_deg_e
    len_deg_n = radius_px / px_per_deg_n

    return calc_compass(image, x, y, len_deg_e, len_deg_n, w)


def calc_compass_center(image, w):
    # calculate center of data
    #x = float(image.shape[1]) / 2.0
    #y = float(image.shape[0]) / 2.0
    #hardcoded values from Jim for trick
    x = 942
    y = 747

    # radius we want the arms to be (approx 1/4 the smallest dimension)
    radius_px = float(min(float(image.shape[1]), float(image.shape[0]))) / 4.0

    return calc_compass_radius(image, x, y, radius_px, w)

def rotate_pt(x_arr, y_arr, theta_deg, xoff=0, yoff=0):
    """
    Rotate an array of points (x_arr, y_arr) by theta_deg offsetted
    from a center point by (xoff, yoff).
    """
    # TODO: use opencv acceleration if available
    a_arr = x_arr - xoff
    b_arr = y_arr - yoff
    cos_t = np.cos(np.radians(theta_deg))
    sin_t = np.sin(np.radians(theta_deg))
    ap = (a_arr * cos_t) - (b_arr * sin_t)
    bp = (a_arr * sin_t) + (b_arr * cos_t)
    return np.asarray((ap + xoff, bp + yoff))

def transform(data_np, flip_x=False, flip_y=False, swap_xy=False):

    # Do transforms as necessary
    if flip_y:
        data_np = np.flipud(data_np)
    if flip_x:
        data_np = np.fliplr(data_np)
    if swap_xy:
        data_np = data_np.swapaxes(0, 1)

    return data_np

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

    if logger is not None:
        logger.debug("rotating with pillow")
    img = PILimage.fromarray(data_np)
    img_rot = img.rotate(theta_deg, resample=False, expand=False,
                         center=(rotctr_x, rotctr_y))
    newdata = np.array(img_rot, dtype=data_np.dtype)
    new_ht, new_wd = newdata.shape[:2]
    assert (wd == new_wd) and (ht == new_ht), \
        Exception("rotated cutout is %dx%d original=%dx%d" % (
            new_wd, new_ht, wd, ht))
    '''
    else:
        if logger is not None:
            logger.debug("rotating with numpy")
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
    '''
    return newdata


def rotate(data_np, theta_deg, rotctr_x=None, rotctr_y=None, pad=20,
           logger=None):

    # If there is no rotation, then we are done
    if math.fmod(theta_deg, 360.0) == 0.0:
        return data_np

    ht, wd = data_np.shape[:2]

    ocx, ocy = wd // 2, ht // 2

    # Make a square with room to rotate
    side = int(math.sqrt(wd**2 + ht**2) + pad)
    new_wd = new_ht = side
    dims = (new_ht, new_wd) + data_np.shape[2:]
    # Find center of new data array
    ncx, ncy = new_wd // 2, new_ht // 2
    '''
        if have_opencl and _use == 'opencl':
            if logger is not None:
                logger.debug("rotating with OpenCL")
            # find offsets of old image in new image
            dx, dy = ncx - ocx, ncy - ocy

            newdata = trcalc_cl.rotate(data_np, theta_deg,
                                       rotctr_x=rotctr_x, rotctr_y=rotctr_y,
                                       clip_val=0, out=None,
                                       out_wd=new_wd, out_ht=new_ht,
                                       out_dx=dx, out_dy=dy)
    '''
    #else:
    # Overlay the old image on the new (blank) image
    ldx, rdx = min(ocx, ncx), min(wd - ocx, ncx)
    bdy, tdy = min(ocy, ncy), min(ht - ocy, ncy)

    # TODO: fill with a different value?
    newdata = np.zeros(dims, dtype=data_np.dtype)
    newdata[ncy - bdy:ncy + tdy, ncx - ldx:ncx + rdx] = \
        data_np[ocy - bdy:ocy + tdy, ocx - ldx:ocx + rdx]

    # Now rotate with clip as usual
    newdata = rotate_clip(newdata, theta_deg,
                          rotctr_x=rotctr_x, rotctr_y=rotctr_y,
                          out=newdata)
    return newdata


def orient(image, w, righthand=False):
        (x, y, xn, yn, xe, ye) = calc_compass_center(image, w)
        degn = math.degrees(math.atan2(xn - x, yn - y))
        # rotate east point also by degn
        xe2, ye2 = rotate_pt(xe, ye, degn, xoff=x, yoff=y)
        dege = math.degrees(math.atan2(xe2 - x, ye2 - y))

        # if right-hand image, flip it to make left hand
        xflip = righthand
        if dege > 0.0:
            xflip = not xflip
        if xflip:
            degn = - degn

        newdata = transform(image, xflip, False, False)
        newdata = rotate(image, degn)
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
    hdu.writeto('rotatedImage.fits')

def main():
    f = fits.open(args.file)
    fitsData = fits.getdata(args.file, ext=0)
    w = wcs.WCS(f[0].header)
    w.wcs.print_contents()
    header = fits.getheader(args.file)
    newdata = orient(fitsData, w)
#    plotIm(newdata)
#    writeFits(header, newdata)


if __name__=="__main__":
    try:
        main()
    except KeyboardInterrupt:
        print('stopping...')