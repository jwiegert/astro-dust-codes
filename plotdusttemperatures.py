# Plot temperature output from thermal simulation. UNDER MODIFICATION!
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


"""
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
"""


#
# Load data (change link to your results folders)
#
modelname = 'testdata/'
#
# Griddistances
#
loaddata = open(modelname+'grid_distances_10degdisc.dat', 'r').readlines()
griddistances = np.array([float(data) for data in loaddata])


sizeofarrays  = open(modelname+'dust_density.inp', 'r').readline()


int(np.loadtxt(modelname+'/dust_density.inp')[1]) #  = np.size(griddistances)
#
# Densities
#
densities     = np.loadtxt(modelname+'/dust_density.inp')[3:]
#
# Temperature
#
temperatures  = np.loadtxt(modelname+'/dust_temperature.dat')[3:]
#
# Extract dusty cells
#
# Downsize arrays by this factor
nres = 100
# Species 1
#
ndust         = np.where(densities[:sizeofarrays] > 0.)[0]
nsmall        = int(np.floor(np.size(ndust)/(1.*nres)))
ndustsmall    = np.zeros(nsmall)
for nn in range(nsmall):
    ndustsmall[nn] = ndust[nres*nn]
ndustsmallA   = ndustsmall.astype(int)
#
# Species 2
#
ndust         = sizeofarrays + np.where(densities[sizeofarrays:] > 0.)[0]
nsmall        = int(np.floor(np.size(ndust)/nres))
ndustsmall    = np.zeros(nsmall)
for nn in range(nsmall):
    ndustsmall[nn] = ndust[nres*nn]
ndustsmallB   = ndustsmall.astype(int)
#
# Standard analytical dust temperature relation ------------------------------ #
#
dustradius              = np.linspace(np.min(inradius)/AUcm,np.max(outradius)/AUcm,100)
dusttemperatureestimate = (0.5 * rstar/AUcm * 1./dustradius)**(2./5.) * startemperature
#                                                                              #
# Plot everything ------------------------------------------------------------ #
#                                                                              #
# plt.savefig('20190215_10to100rstar_spiral_temperature_manyphotonsinfinegrid.png')
#
plt.ion()
#
plt.figure('Temperature of dust cells', figsize=(6, 6))
rc('xtick.major',size=8)
rc('xtick.minor',size=4)
rc('ytick.major',size=8)
rc('ytick.minor',size=4)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#
# Species 1 (sphere)
#
plt.subplot(1,2,1)
plt.plot([np.min(inradius)/AUcm,np.max(outradius)/AUcm],[Tin ,Tin] ,'k')
plt.plot([np.min(inradius)/AUcm,np.max(outradius)/AUcm],[Tout,Tout],'k')
plt.plot(griddistances[ndustsmallA]/AUcm,temperatures[ndustsmallA],'b.',markersize=2)
#
# And analytical estimate
#
plt.plot(dustradius,dusttemperatureestimate,'k--')
#
plt.annotate(r'Species1: sphere',xy=(80,300))
#plt.xlim(0,np.max(griddistances[ndustsmallA])/AUcm+10)
#plt.ylim(0,np.max(temperatures[ndustsmallA])+10)
plt.xlabel(r'Radius (AU)',fontsize=18)
plt.ylabel(r'Temperature (K)',fontsize=18)
#plt.xscale('log')
#plt.yscale('log')
#
# Species 2 (disc)
#
plt.subplot(1,2,2)
plt.plot([np.min(inradius)/AUcm,np.max(outradius)/AUcm],[Tin ,Tin] ,'k')
plt.plot([np.min(inradius)/AUcm,np.max(outradius)/AUcm],[Tout,Tout],'k')
plt.plot(griddistances[ndustsmallB-sizeofarrays]/AUcm,temperatures[ndustsmallB], 'r.',markersize=2)
#
# And analytical estimate
#
plt.plot(dustradius,dusttemperatureestimate,'k--')
#
plt.annotate(r'Species2: disc',xy=(80,300))
#plt.xlim(0,np.max(griddistances[ndustsmallB-sizeofarrays])/AUcm+10)
#plt.ylim(0,np.max(temperatures[ndustsmallB])+10)
plt.xlabel(r'Radius (AU)',fontsize=18)
plt.ylabel(r'Temperature (K)',fontsize=18)
#plt.xscale('log')
#plt.yscale('log')
#
#
#

