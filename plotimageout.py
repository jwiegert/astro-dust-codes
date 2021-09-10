# Loads and plots the radmc3d-computed SEDs with and without star.
# ============================================================================ #
#
import numpy as np
import scipy as s
from scipy import ndimage
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['serif']})
rc('text', usetex=True)
#
# Useful Constants
# 
from constants import *
from inputdata import *
#
# useful function(s)
#
import settickslabels
#
# ============================================================================ #
#
# Load the data. Flux density is in erg cm-2 s-1 Hz-1 ster-1.
# Yes, manually change the filename :P
# This variable is since I often have many images and want to switch between
# them easily.
#
filename   = 'image50um_edgeon.out'
#
imagedata  = open(filename, 'r').readlines()
npix       = int(imagedata[1].split()[0])
pixsizeau  = float(imagedata[3].split()[0]) / AUcm # Pixel size in AU
pixsize    = pixsizeau / distance           # and in asec
wavelength = float(imagedata[4])
image1d    = np.zeros(npix*npix)
for n,data in enumerate(imagedata[6:-2]):
    image1d[n] = float(data)
#
# Set up image axis
#
sizeau    = npix   * pixsizeau # Length of side of image in AU
sizemas   = sizeau / distance  # Length of side of image in mas
axisplot  = [0.5*sizeau,-0.5*sizeau,-0.5*sizeau,0.5*sizeau]
axisrange = np.linspace(-0.5*sizeau,0.5*sizeau,npix)
#
# Recalculate flux to Jy asec-2  at distance to the star
# 1sec2 = 2.35044305391e-11 ster
#
image1d    = image1d * 1.e23 * asecsqster * pixsize**2.
#
# Extract limits
#
fluxmin,fluxmax = np.min(image1d),np.max(image1d)
logfluxrange    = np.logspace(np.log10(fluxmax)-10.,np.log10(fluxmax),10)
#
# Extract 2D image from input image data
#
image2d    = np.zeros((npix,npix))
logimage2d = np.zeros((npix,npix))
npixtot    = int(npix*npix)
nx,ny      = 0,0
for nn in range(npixtot):
    image2d[nx,ny] = image1d[nn]
    if image1d[nn] <= logfluxrange[0]:
        logimage2d[nx,ny] = np.log10(logfluxrange[0])
    else:
        logimage2d[nx,ny] = np.log10(image1d[nn])
    nx += 1
    if nx == npix:
        nx = 0
        ny += 1
#
# Apply filter to emulate telescope resolution
# Sigma in pixels (FWHM is 5.6 asec, FWHM ~ 2 sqrt 2 ln sigma )
#
resolution     = 5.6 / 2.355
filtimage2d    = s.ndimage.gaussian_filter(image2d, sigma=resolution)
filtimagemax   = np.max(filtimage2d)
logfiltimage2d = s.ndimage.gaussian_filter(logimage2d, sigma=resolution)
logmax         = np.log10(filtimagemax)
#
# Plot and print outputs. ==================================================== #
#
plt.ion()
#
# Image
#
plt.figure('Simulated original images')
for nn in range(6):
    plt.subplot(2,3,1+nn)
    plt.imshow(image2d, origin='lower', interpolation='nearest', vmin=fluxmin, vmax=np.logspace(-6,0,6)[nn]*fluxmax, extent=axisplot)
    cb = plt.colorbar(orientation = 'vertical',shrink=0.6,pad=0.15)
    cb.set_label(label = r'Flux density, (Jy/arcsec$^2$)',fontsize= 18)
    cb.ax.tick_params(labelsize=18)
    settickslabels.settickslabels(xlabel="Offset (AU)", ylabel="Offset (AU)", xscale="lin", yscale="lin")
#
# Filtered image
#
plt.figure('Simulated convolved image')
plt.imshow(filtimage2d, origin='lower', interpolation='nearest', vmin=fluxmin, vmax=filtimagemax, extent=axisplot)
cb = plt.colorbar(orientation = 'vertical',shrink=0.6,pad=0.15)
cb.set_label(label = r'Flux density, (Jy/arcsec$^2$)',fontsize= 18)
cb.ax.tick_params(labelsize=18)
settickslabels.settickslabels(xlabel="Offset (AU)", ylabel="Offset (AU)", xscale="lin", yscale="lin")
#
# Logged Filtered image
#
plt.figure('Simulated convolved logged image')
plt.imshow(logfiltimage2d, origin='lower', interpolation='nearest', vmin=np.log10(logfluxrange[0]), vmax=logmax, extent=axisplot)
cb = plt.colorbar(orientation = 'vertical',shrink=0.6,pad=0.15)
cb.set_label(label = r'Log10 of flux density, (Jy/arcsec$^2$)',fontsize= 18)
cb.ax.tick_params(labelsize=18)
settickslabels.settickslabels(xlabel="Offset (AU)", ylabel="Offset (AU)", xscale="lin", yscale="lin")
#
# Filtered image contours
#
plt.figure('Simulated convolved contours')
plt.contour(axisrange,axisrange,filtimage2d,logfluxrange)
settickslabels.settickslabels(xlabel="Offset (AU)", ylabel="Offset (AU)", xscale="lin", yscale="lin")
