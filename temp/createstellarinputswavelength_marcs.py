# Loads model data, observed SED and spectra and transforms units to Jy.       #
# Creates input data for R3D from the stellar MARCS-model, wavelength grid,    #
# and saves files used for plotting of results from R3D.                       #
#                                                                              #
#==============================================================================#
import sys                                                                     #
sys.path.append('/home/joachim/programexempel/python/pythonstartup')           #
from pythonstartup import *                                                    #
#==============================================================================#
from inputdata import *
#
# Load data. ----------------------------------------------------------------- #
#
# Observed stellar SED (0: wavelength, 1: Flux jy, 2: error)
#
starsed = np.loadtxt('martin/iras17521_sed_jansky.dat')
#
# PACS-spectrum (0: wavelength, 1: Flux Jy, 2: error Jy)
#
starpacsspectrum = np.loadtxt('iras17521_69_pacsspectrum.dat')
#
# ISO(?) spectrum (0: wavelength, 1: Flux Jy, 2: error Jy)
#
starirsspectrum = np.loadtxt('iras17553_irsspectrum.dat')
#
# Create stars.inp ----------------------------------------------------------- #
#
# MARCS Starmodel, extracted by me with the same properties as Martin's star.
# Already in Jy.
#
starmodel = np.zeros((nlambda,2))
starmodel[:,1]   = np.loadtxt('marcs/stars.dat')
starmodel[:,0]   = np.loadtxt('marcs/wavelength.dat')
#
# Compute luminosity of the model SED and a lum-corr-factor to correct for the
# lack of non-point-source-handling of r3d.
#
# Uses a trapezoidal integral.
#
sedintegral         = 0.
for nn in range(nlambda - 1):
    sedintegral     = sedintegral + \
                  0.5*(starmodel[nn,1] + starmodel[nn+1,1])*1e-26 * \
                      (c/(starmodel[nn,0]*1e-6) - c/(starmodel[nn+1,0]*1e-6))
lumcorr             = luminosity / (4.*np.pi*(distance*pc)**2. * sedintegral)
lumdistanceunitcorr = 1.e-23 * distance**2. * lumcorr
print('Unit, distance, and luminosity correction factor (from Jy at observed distance to spectra as used by R3d):', lumdistanceunitcorr)
#
# 1. Convert flux units to R3D's units (in erg/cm^2/s/Hz at 1pc)
#
#       1 Jy = 1e-23 erg/cm^2/s/Hz
#
# Also multiply by 4 pi rstar**2, because I want the flux density to be at the star's surface. rstar is in cm, and distance is in pc > rstar/(pc*100)
#
# Drop the 4pi because the distance correction factor is 4pi*distance(pc)^2 / 4pi*1pc^2
#
starr3d = starmodel[:,1] * lumdistanceunitcorr
#
# 2. Write stars.inp and don't forget to include wavelength grid in it.
#
f = open('stars.inp','w')
f.writelines(['2\n',\
              str(nstar).strip('[]'),' ',str(np.size(starmodel[:,1])).strip('[]'),'\n',\
              str(rstar).strip('[]'),' ',str(mstar).strip('[]'),' ',str(coordinatesstars[0]).strip('[]'),' ',str(coordinatesstars[0]).strip('[]'),' ',str(coordinatesstars[0]).strip('[]'),'\n'])
for nn in range(np.size(starmodel[:,1])):
    f.writelines([str(starmodel[nn,0]).strip('[]'),'\n'])
for nn in range(np.size(starmodel[:,1])):
    f.writelines([str(starr3d[nn]).strip('[]'),'\n'])
f.close()
#
# 3. Write wavelength grid (used for plotting)
#
np.savetxt('wavelength.dat',starmodel[:,0],fmt='%f')
#
# Create R3D's wavelength file, wavelength_micron.inp
#
f = open('wavelength_micron.inp','w')
f.writelines([str(np.size(starmodel[:,0])).strip('[]'),'\n'])
for nn in range(np.size(starmodel[:,0])):
    f.writelines([str(starmodel[nn,0]).strip('[]'),'\n'])
f.close()
#
# Plot data to check that everything is correct. ----------------------------- #
plt.ion()
#
# Plot MARCS model
#
plt.plot(starmodel[:,0],starmodel[:,1])
#
# Plot observed SED
#
plt.errorbar(starsed[:,0], starsed[:,1], yerr=starsed[:,2],\
             fmt='o', markeredgecolor='k', color='k')
#
# Plot spectra
#
plt.plot(starpacsspectrum[:,0],starpacsspectrum[:,1],'g--')
plt.plot(starirsspectrum[:,0] ,starirsspectrum[:,1] ,'y--')
#
# TEMP: compare with a black body of 2600K
#
testsed          = 1e23/(distance + rstar/(pc*100))**2. * \
                   np.loadtxt('r3dresults/spectrumincl00_iras17521_blackbody.out')[:,1]
testwavelength   = np.loadtxt('r3dresults/spectrumincl00_iras17521_blackbody.out')[:,0]
plt.plot(testwavelength,testsed,'y--')
#
# Axes properties
#
plt.xscale('log')
plt.xlim(0.1,2000.)
plt.xticks(fontsize=18)
plt.xlabel('Wavelength, $\lambda$ ($\mu$m)',fontsize=18)
plt.yscale('log')
plt.ylim(1e-12,1.e2)
plt.yticks(fontsize=18)
plt.ylabel(r'Flux density, $F_\nu$ (Jy)',fontsize=18)


