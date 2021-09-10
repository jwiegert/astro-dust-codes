# Load and change total dust mass of dust_density
#
import sys                                                                     #
sys.path.append('/home/joachim/programexempel/python/pythonstartup')           #
from pythonstartup import *                                                    #
#==============================================================================#
#from inputdata import *
#
# Load original density file
#
temp     = np.loadtxt('dust_density.inp')
density1 = temp[3:]
nleafs   = int(temp[1])
nrspec   = int(temp[2])
del temp
print('')
print('Loaded data.')
#
# Extract original and wanted mass (run r3d to see current mass)
#
mass1 = 1e-8                    # Current mass
mass2 = 2e-8                    # Wanted mass
#
# Modify mass
#
density2 = density1 * mass2/mass1
#
# Save new file
#
print('Current mass:',mass1)
print('Saving new mass:',mass2)
f = open('dust_density.inp', 'w')
f.writelines(['1\n',\
              str(nleafs),'\n',\
              str(nrspec),'\n'])
for nspecies in range(nrspec):
    for nn in range(nleafs):
        f.writelines([str(density2[nspecies*nleafs + nn]).strip('[]'),'\n'])
f.close()
print('Done')
#
# Make a noise also :)
#
os.system('play --no-show-progress --null --channels 1 synth %s sine %f' %( 0.1, 400))

