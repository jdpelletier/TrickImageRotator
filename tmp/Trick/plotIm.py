import argparse
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
from matplotlib import colors


parser = argparse.ArgumentParser(description="Plot fits image ",
                         usage="plotIm.py file")

parser.add_argument("file", help="file to load")


args = parser.parse_args()

plt.style.use(astropy_mpl_style)


image_data = fits.getdata(args.file, ext=0)

plt.figure()
plt.imshow(image_data, cmap='gray', norm=colors.LogNorm())
plt.colorbar()
plt.show()
