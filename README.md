# TrickImageRotator, UNDER DEVELOPMENT
Tool for OA's/observers to have a rotated full frame processed Trick image.

## Prerequisites for running
Only know it runs on k1aoserver-new right now, will update as this changes

## Usage
### For ginga (primary) version:
FIRST you must set the LD_LIBRARY_PATH to get ginga to work.  Run this on k1aoserver-new:

setenv LD_LIBRARY_PATH /usr/local/anaconda/lib:/kroot/rel/default/lib:/usr/lib:/usr/lib64

This will set the variable for your session, and must be ran each time you start
a new session.  After this, you can run ginga and use this tool:

On k1aoserver-new, in the /home/k1obsao/jpelletier directory, run the command 
"kpython3 gignaTrickDisplay.py". This will open the gui, where you can open images
manually or start the scanning function, which will search the current nightly 
directory for images as the come in.  When opening manually, the window will
also be set to the nightly directory so you won't have to search around for images.

### For ds9 (backup) version:
On k1aoserver-new, run the command "kpython3 imageDisplay.py". The script will start 
scanning the directory /net/k1aoserver/k1aodata/nightly for new images. Once a new
image is found, the data is pulled and processed, then displayed with ds9_80.

## Current functions
- Rotates viewer to get the image north up, east left (alledgedly....)
  - Because the viewer is rotating, pixel information is retained
- Displays current ROI spot
- Adds compass rose and scale for user (scale only in ds9 right now)
- Ginga should be the primary display, but ds9 can work as a backup
  - ds9 rotates the image, so pixel information is not retained

## Future additions
-Boxes showing where OSIMG and OSPEC are looking
