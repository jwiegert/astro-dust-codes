# Creates mass absorption and scattering input files for R3D from Martin's     #
# data on the opacity Qabs and Qscat.                                          #
#                                                                              #
#==============================================================================#
import sys                                                                     #
sys.path.append('/home/joachim/programexempel/python/pythonstartup')           #
from pythonstartup import *                                                    #
#==============================================================================#
from inputdata import *
#
# Load Qabs and Qscat
#
wavelengthFost   = np.loadtxt('martin/dhs0.7_a0.50ForsteriteSuto200.dat')[:,0]
QabsFost         = np.loadtxt('martin/dhs0.7_a0.50ForsteriteSuto200.dat')[:,1]
QscatFost        = np.loadtxt('martin/dhs0.7_a0.50ForsteriteSuto200.dat')[:,2]
wavelengthAmIron = np.loadtxt('martin/dhs0.7_a0.50MgFeSiO4100AlOx0cdeFe10.dat')[:,0]
QabsAmIron       = np.loadtxt('martin/dhs0.7_a0.50MgFeSiO4100AlOx0cdeFe10.dat')[:,1]
QscatAmIron      = np.loadtxt('martin/dhs0.7_a0.50MgFeSiO4100AlOx0cdeFe10.dat')[:,2]
#
# Calculate mass extinction, kappa, in cm2/g
#
# KappaExt = Qext * 3/(4*a*rho)
#   a = grain radius in cm
# rho = grain density in g/cm3
#
# (a is in um here:)
# Fosterite (Mg2SiO4): 
# a=  0.500 CrOlivSuto200  rho= 2.11
#
# Amorphous+Iron:
# MgFeSiO4 100 AlO x0cde Fe10.dat
#a=  0.500 MgFeSiO4      rho= 2.14        
#a=  0.500 AlOx-compac   rho= 2.60        NOT INCLUDED
#a=  0.200  CrOlivine    rho= 2.11        NOT INCLUDED
#a=  0.010 IronPollack   rho= 7.8                     
# massfractions= 0.91 0.00 0.00 0.09 > rho= 2.65
# 
agrain          = 0.5 * 1e-4     # Grain radius in cm, 0.5um
rhoFost         = 2.11           # Grain density in g/cm3
rhoAmIron       = 2.65           # Grain density in g/cm3
corrFost        = 3./(4.*agrain*rhoFost)
corrAmIron      = 3./(4.*agrain*rhoAmIron)
#
kappaAbsFost    = QabsFost    * corrFost
kappaScatFost   = QscatFost   * corrFost
kappaAbsAmIron  = QabsAmIron  * corrAmIron
kappaScatAmIron = QscatAmIron * corrAmIron
#
# Mean scattering angle ⟨cos(θ)⟩ = g(lambda), not used, set to 0.
#
scatangleFost   = np.zeros(np.size(wavelengthFost))
scatangleAmIron = np.zeros(np.size(wavelengthAmIron))
#
# Check the files for duplicates
#
for nn in range(np.size(wavelengthFost)-1):
    if wavelengthFost[nn] == wavelengthFost[nn+1]:
        print('Fosterite: ', wavelengthFost[nn])
for nn in range(np.size(wavelengthAmIron)-1):
    if wavelengthAmIron[nn] == wavelengthAmIron[nn+1]:
        print('amIron: ', wavelengthAmIron[nn])

#
# Write files for R3D -------------------------------------------------------- #
#
# 1. dustopac file that lists dust species.
#
f = open('dustopac.inp','w')
f.writelines(['2\n',\
              str(nrspec),'\n',\
              '-----------------------------\n',\
              '1\n',\
              '0\n',\
              'amIron\n',\
              '-----------------------------\n',\
              '1\n',\
              '0\n',\
              'Fosterite\n',\
              '-----------------------------'])
f.close()
#
# 2. Fosterite file
#
f = open('dustkappa_Fosterite.inp','w')
f.writelines(['# Fosterite (Mg2SiO4): (Ask Martin Groenewegen for the reference)\n',\
              '# a=0.500 um, CrOlivSuto200  rho= 2.11 g/cm3\n',\
              '#\n',\
              '# Columns are:\n',\
              '# Wavelength(um), Kappa_Abs(cm2/g), Kappa_Scat(cm2/g), g(lambda) = mean scattering angle\n',\
              '2\n',\
              str(np.size(wavelengthFost)),'\n'])
for nn in range(np.size(wavelengthFost)):
    f.writelines([str(wavelengthFost[nn]).strip('[]'),'    ',\
                  str(kappaAbsFost[nn]).strip('[]'),'    ',  \
                  str(kappaScatFost[nn]).strip('[]'),'    ', \
                  str(scatangleFost[nn]).strip('[]'),'\n'])
f.close()
#
# 3. Amorphous Iron file
#
f = open('dustkappa_amIron.inp','w')
f.writelines(['# Amorphous+Iron: (Ask Martin Groenewegen for the reference)\n',\
              '# MgFeSiO4 100 AlO x0cde Fe10.dat\n',\
              '#a=  0.500um MgFeSiO4      rho= 2.14g/cm3        \n',\
              '#a=  0.500um AlOx-compac   rho= 2.60g/cm3        NOT INCLUDED\n',\
              '#a=  0.200um  CrOlivine    rho= 2.11g/cm3        NOT INCLUDED\n',\
              '#a=  0.010um IronPollack   rho= 7.8 g/cm3                    \n',\
              '# massfractions= 0.91 0.00 0.00 0.09 > rho= 2.65g/cm3\n',\
              '#\n',\
              '# Columns are:\n',\
              '# Wavelength(um), Kappa_Abs(cm2/g), Kappa_Scat(cm2/g), g(lambda) = mean scattering angle\n',\
              '2\n',\
              str(np.size(wavelengthAmIron)),'\n'])
for nn in range(np.size(wavelengthAmIron)):
    f.writelines([str(wavelengthAmIron[nn]).strip('[]'),'    ',\
                  str(kappaAbsAmIron[nn]).strip('[]'),'    ',  \
                  str(kappaScatAmIron[nn]).strip('[]'),'    ', \
                  str(scatangleAmIron[nn]).strip('[]'),'\n'])
f.close()
#
# Plot to check -------------------------------------------------------------- #
#
plt.ion()
plt.subplot(1,2,1)
plt.plot(wavelengthFost,QabsFost,'b')
plt.plot(wavelengthFost,QscatFost,'b--')
plt.plot(wavelengthAmIron,QabsAmIron,'r')
plt.plot(wavelengthAmIron,QscatAmIron,'r--')
#
plt.xscale('log')
plt.xlim(0.1,2000.)
plt.xticks(fontsize=18)
plt.xlabel('Wavelength, $\lambda$ ($\mu$m)',fontsize=18)
plt.yscale('log')
plt.yticks(fontsize=18)
#
#
#
plt.subplot(1,2,2)
plt.plot(wavelengthFost,kappaAbsFost,'b')
plt.plot(wavelengthFost,kappaScatFost,'b--')
plt.plot(wavelengthAmIron,kappaAbsAmIron,'r')
plt.plot(wavelengthAmIron,kappaScatAmIron,'r--')
#
plt.xscale('log')
plt.xlim(0.1,2000.)
plt.xticks(fontsize=18)
plt.xlabel('Wavelength, $\lambda$ ($\mu$m)',fontsize=18)
plt.yscale('log')
plt.yticks(fontsize=18)

oldkappa = np.loadtxt('kappaext_mg2sio4-10nm-compare.dat')
plt.subplot(1,2,2)
plt.plot(oldkappa[:,0],oldkappa[:,1],'y')
plt.plot(oldkappa[:,0],oldkappa[:,2],'y--')




