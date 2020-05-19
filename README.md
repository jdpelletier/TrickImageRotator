# TrickImageRotator, UNDER DEVELOPMENT
Tool for OA's/observers to have a rotated full frame processed Trick image.

# Prerequisites for running
Only know it runs on k1aoserver-new right now, will update as this changes

# Usage
On k1aoserver-new, run the command "kpython3 imageDisplay.py". The script will start 
scanning the directory /net/k1aoserver/k1aodata/nightly for new images. Once a new
image is found, the data is pulled and processed, then displayed with ds9_80.

# Current functions
>Rotates image to north up, east left (alledgedly....)
>Displays current ROI spot
>Adds compass rose and scale for user

# Future additions
>Boxes showing where OSIMG and OSPEC are looking
>Ditching ds9 for python built gui
