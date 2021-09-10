# Extracts down-sampled versions of the observed spectra for chi2-comparisons.
#
#==============================================================================#
# Header for python scripts (remove interpolate if you don't need it ==========#
import sys                                                                     #
sys.path.append('/home/joachim/programexempel/python/pythonstartup')           #
from pythonstartup import *                                                    #
#==============================================================================#
#
from inputdata import *
#
# Load wavelength grid
#
wavelength     = np.loadtxt('wavelength.dat')
wavelengthsize = np.size(wavelength)
#
# Sample these wavelengths from the spectras and save as sampled spectras
# 1. sws spectrum
# [:,0] = wavelength (um)
# [:,1] = Flux (Jy)
# [:,2] = error (Jy)
# [:,3] = error (Jy) incl normalization error
#
# 1. PACS-spectrum (0: wavelength, 1: Flux Jy, 2: error Jy)
#
pacsspectrum     = np.loadtxt('iras17521_69_pacsspectrum.dat')
wavelengthpacs   = pacsspectrum[:,0]
#
# Find the number of r3dwavelengths within the spectra
#
pacslength       = np.size(np.where(wavelength > wavelengthpacs[0])[0]) - \
                   np.size(np.where(wavelength > wavelengthpacs[-1])[0])
pacsindex        = np.zeros(pacslength)
#
# 2. ISO(?) IRAS(?) IRS(?) spectrum (0: wavelength, 1: Flux Jy, 2: error Jy)
#
irsspectrum = np.loadtxt('iras17553_irsspectrum.dat')
wavelengthirs    = irsspectrum[:,0]
#
# Find the number of r3dwavelengths within the spectra
#
irslength        = np.size(np.where(wavelength > wavelengthirs[0])[0]) - \
                   np.size(np.where(wavelength > wavelengthirs[-1])[0])
irsindex         = np.zeros(irslength)
#
# Extract the indeces of the corresponding wavelengths from the spectra
#
counterpacs      = 0
counterirs       = 0
for nn in range(wavelengthsize):
    if      wavelength[nn] > wavelengthpacs[ 0]\
    and     wavelength[nn] < wavelengthpacs[-1]:
        pacsindex[counterpacs] = np.argwhere(wavelengthpacs > wavelength[nn])[0][0]
        counterpacs        = counterpacs + 1
    if      wavelength[nn] > wavelengthirs[ 0]\
    and     wavelength[nn] < wavelengthirs[-1]:
        irsindex[counterirs]  = np.argwhere(wavelengthirs > wavelength[nn])[0][0]
        counterirs         = counterirs + 1
#
# Use the indeces to extract and save the sampled data (wavelength, data, error)
#
pacsindex       = pacsindex.astype(int)
pacssample      = np.zeros((pacslength,3))
pacssample[:,0] = pacsspectrum[pacsindex,0]
pacssample[:,1] = pacsspectrum[pacsindex,1]
pacssample[:,2] = pacsspectrum[pacsindex,2]
np.savetxt('iras17521_69_pacsspectrum_sampled.dat',pacssample,fmt='%f')
irsindex        = irsindex.astype(int)
irssample       = np.zeros((irslength,3))
irssample[:,0]  = irsspectrum[irsindex,0]
irssample[:,1]  = irsspectrum[irsindex,1]
irssample[:,2]  = irsspectrum[irsindex,2]
np.savetxt('iras17553_irsspectrum_sampled.dat',irssample,fmt='%f')
#    
# Plot to check -------------------------------------------------------------- #
#  sws spectrum
#
plt.ion()
plt.figure(r'IRS and PACS')
plt.plot(wavelengthirs,irsspectrum[:,1],'orange',linewidth=2)
plt.plot(irssample[:,0],irssample[:,1],'m')
plt.plot(irssample[:,0],irssample[:,1]+irssample[:,2],'m:')
plt.plot(irssample[:,0],irssample[:,1]-irssample[:,2],'m:')
plt.xscale('log')
plt.yscale('log')
#
# Pacs spectrum
#
plt.plot(wavelengthpacs,pacsspectrum[:,1],'cyan',linewidth=2)
plt.plot(wavelengthpacs,pacsspectrum[:,1],'cyan',linewidth=2)
plt.plot(pacssample[:,0],pacssample[:,1],'b')
plt.plot(pacssample[:,0],pacssample[:,1]+pacssample[:,2],'b:')
plt.plot(pacssample[:,0],pacssample[:,1]-pacssample[:,2],'b:')

