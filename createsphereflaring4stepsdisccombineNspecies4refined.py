# This file creates a SPHERE and a DISC with a fraction of the sphere's mass
# and the GRID for the combination (same grid as for sphere and disc but
# combined). It assumes N dust species and FOUR grid refinements.
#
# The disc here is flaring.
# With a flaring relation (height as function of radius):
#      h(r) = flaringconstant * r
#
# The flaring is "square'ey" so, at each
#      R = 2*gridsmaller/0.1 = 20.*gridsmaller
# The disc doubles in height.


# ============================================================================ #
# Import packages
import os
import numpy as np
import scipy as s
from scipy.integrate import quad
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['serif']})
rc('text', usetex=True)


# Import constants and properties.
from constants import *
from inputdata import *

# ============================================================================ #
# Additional function(s)
def movecoordinates(nx,ny,nz):
    nx     += 1
    if nx  == nxyz:
        nx =  0
        ny += 1
    if ny  == nxyz:
        nx =  0
        ny =  0
        nz += 1
    return nx,ny,nz

# ============================================================================ #

savegrid = input("Do you want to save the grid distances file? y/n: ")
plotgrid = input("Do you want to plot the grid? y/n: ")


# List of inputs
# 
# sphereratio
# denspower (list with 2 elements, sphere and disc, radial power) defaults; [-2,-2]
# totaldustmass
# flareangle 
#        (flaring in degrees, gives h(r) = discflare * r)
#        (flare angle is angle from mid-plane to disc edge)
#        (tan(flareangle) = 2*discflare )
#        (discflare = 2. * np.tan(0.017453292519943295 * flareangle) )
#        (0.0174... = np.pi/180)
#        (h(r) = r * 2 * tan(flareangle) = r * discflare)


# inradius (list with 2 elements, sphere and disc)
# outradius (list with 2 elements, sphere and disc)
# save grid distances? default Yes
# plot the grid? default No


# discmasscorr = beror på discradie och flaring




#def create_envelope():



#
# Reminder concerning disc height:
# h(r) = discwidth(r)*2
# h(r) = discflare*r
# -> r = 2*discwidth(r)/discflare
#
# Sphere and disc properties
#
discratio      = (1. - sphereratio) * discmasscorr
#
# Smallest cube will always have its edges in the grid origin. It must be a
# at least a quarter of the inner hole as the largest. Otherwise the disc will
# be missed and hole will be too "squarey".
#
smallcubesize4 = 0.2 * np.min(inradius)
#
print(f"""
All outputs are also saved in grid_info.txt

                 Sphere between: {inradius[0]/AUcm} AU to {outradius[0]/AUcm} AU
                 Sphere between: {inradius[0]/AUcm} AU to {outradius[0]/AUcm} AU
                   Disc between: {inradius[1]/AUcm} AU to {outradius[1]/AUcm} AU
  Approx. disc flaring relation: {discflare} * r
              Sphere mass ratio: {sphereratio}
                  Input mass is: {totaldustmass/1.989e33/spheremasscorr} Msol
Sphere, species1 density distr.: propto r^{denspower[0]}
  Disc, species2 density distr.: propto r^{denspower[1]}""")
#                                                                              #
# Compute normalisation densities. =========================================== #
#                                                                              #
spherezerodensity = 0.
disczerodensity   = 0.
for nspecies in range(nrspec):
    intR = s.integrate.quad(lambda x: \
        x**(2.+denspower[nspecies]), inradius[nspecies], outradius[nspecies])
    if nspecies == 0:
        spherezerodensity = \
            sphereratio * totaldustmass * \
            inradius[nspecies]**denspower[nspecies] / (4.*np.pi * intR[0])
    if nspecies == 1:
        disczerodensity = \
            discratio * totaldustmass * \
            inradius[nspecies]**denspower[nspecies] / (4.*np.pi * intR[0])
totalzerodensity = spherezerodensity + disczerodensity
print(f"""
Density percentage
sphere/(disc+sphere) = {spherezerodensity / (disczerodensity + spherezerodensity)}
disc  /(disc+sphere) = {disczerodensity / (disczerodensity + spherezerodensity)}""")
#
# Grid properties. =========================================================== #
#
# So... the first "basecubesize", and "nxyz" is 
# basecubesize = 2**4 * smallcubesize4
#
nrefines       = 4
basecubesize   = 16. * smallcubesize4
nxyz           = int(np.ceil(2.*np.max(outradius) / basecubesize))
#
# makes sure that nxyz is even, so that the disc is centrered at z=0.
#
if nxyz%2 != 0:
    nxyz += 1
nbasecubes     = int(nxyz * nxyz * nxyz)
gridedge       = nxyz * basecubesize
#
# Define basegrid and a distance to center array
#
gridcourners   = np.linspace(-gridedge*0.5,gridedge*0.5,nxyz+1)
griddist       = np.zeros(nxyz)
for nn in range(nxyz):
    griddist[nn]  = 0.5 * (gridcourners[nn] + gridcourners[nn+1])
#
# Re-compute the grid properties from this
#
basecubesize   = gridcourners[1] - gridcourners[0]
smallcubesize1 = 0.5 * basecubesize
smallcubesize2 = 0.5 * smallcubesize1
smallcubesize3 = 0.5 * smallcubesize2
smallcubesize4 = 0.5 * smallcubesize3
gridcorrx1     = np.array([-1,1,-1,1,-1,1,-1,1]) * 0.5 * smallcubesize1
gridcorry1     = np.array([-1,-1,1,1,-1,-1,1,1]) * 0.5 * smallcubesize1
gridcorrz1     = np.array([-1,-1,-1,-1,1,1,1,1]) * 0.5 * smallcubesize1
gridcorrx2     = np.array([-1,1,-1,1,-1,1,-1,1]) * 0.5 * smallcubesize2
gridcorry2     = np.array([-1,-1,1,1,-1,-1,1,1]) * 0.5 * smallcubesize2
gridcorrz2     = np.array([-1,-1,-1,-1,1,1,1,1]) * 0.5 * smallcubesize2
gridcorrx3     = np.array([-1,1,-1,1,-1,1,-1,1]) * 0.5 * smallcubesize3
gridcorry3     = np.array([-1,-1,1,1,-1,-1,1,1]) * 0.5 * smallcubesize3
gridcorrz3     = np.array([-1,-1,-1,-1,1,1,1,1]) * 0.5 * smallcubesize3
gridcorrx4     = np.array([-1,1,-1,1,-1,1,-1,1]) * 0.5 * smallcubesize4
gridcorry4     = np.array([-1,-1,1,1,-1,-1,1,1]) * 0.5 * smallcubesize4
gridcorrz4     = np.array([-1,-1,-1,-1,1,1,1,1]) * 0.5 * smallcubesize4
#
# Compute disc thicknesses
# 0.75 is from 0.5*1.5
# 0.25 is from 0.5*0.5
# 0.5 because vertical limit is half thickness, 1.5 because I want the height
# at the middle of the final radial range, 0.5 because I want the middle of the
# other ranges.
#
radiallimit1   = 0.5   * np.max(outradius)
radiallimit2   = 0.25  * np.max(outradius)
radiallimit3   = 0.125 * np.max(outradius)
radiallimit4   = 0.0625* np.max(outradius)
#
verticallimit0 = 0.75*discflare*radiallimit1
verticallimit1 = 0.25*discflare*radiallimit1
verticallimit2 = 0.25*discflare*radiallimit2
verticallimit3 = 0.25*discflare*radiallimit3
verticallimit4 = 0.25*discflare*radiallimit4
#
# Compute number of cells inside the vertical limits for the real vertical limit
#
numbcells0 = np.ceil(verticallimit0 /  basecubesize )
numbcells1 = np.ceil(verticallimit1 / smallcubesize1)
numbcells2 = np.ceil(verticallimit2 / smallcubesize2)
numbcells3 = np.ceil(verticallimit3 / smallcubesize3)
numbcells4 = np.ceil(verticallimit4 / smallcubesize4)
#
# Compute a better vertical limit in whole cells
#
verticallimit0 = numbcells0 *  basecubesize
verticallimit1 = numbcells1 * smallcubesize1
verticallimit2 = numbcells2 * smallcubesize2
verticallimit3 = numbcells3 * smallcubesize3
verticallimit4 = numbcells4 * smallcubesize4
#
# Smallest size allowed is twice the grid cell size.
# These ifs makes sure that this happens.
#
if  verticallimit0 <  basecubesize :
    verticallimit0 =  basecubesize
if  verticallimit1 < smallcubesize1 :
    verticallimit1 = smallcubesize1
if  verticallimit2 < smallcubesize2 :
    verticallimit2 = smallcubesize2
if  verticallimit3 < smallcubesize3 :
    verticallimit3 = smallcubesize3
if  verticallimit4 < smallcubesize4 :
    verticallimit4 = smallcubesize4
#
# More output
#
print(f"""
     Final smallest grid size: {smallcubesize4/AUcm} AU
      Third refined grid size: {smallcubesize3/AUcm} AU
     Second refined grid size: {smallcubesize2/AUcm} AU
      First refined grid size: {smallcubesize1/AUcm} AU
            Largest grid size: {basecubesize/AUcm} AU

Disc half thicknesses:
 Inside Radius = {radiallimit4/AUcm} Thickness: {verticallimit4/AUcm} AU, Number of 4th refined cells {numbcells4}
 Inside Radius = {radiallimit3/AUcm} Thickness: {verticallimit3/AUcm} AU, Number of 3rd refined cells {numbcells3}
 Inside Radius = {radiallimit2/AUcm} Thickness: {verticallimit2/AUcm} AU, Number of 2nd refined cells {numbcells2}
 Inside Radius = {radiallimit1/AUcm} Thickness: {verticallimit1/AUcm} AU, Number of 1st refined cells {numbcells1}
Outside Radius = {radiallimit1/AUcm} Thickness: {verticallimit0/AUcm} AU, Number of base cells        {numbcells0}
""")
#
# Define refinment arrays. =================================================== #
# At each level we check first the sphere and then the disc and add refinements
# where it's appropiate.
#
# First refinements
#
refinearray1      = np.zeros(nbasecubes)
nx,ny,nz          = 0,0,0
for nn in range(nbasecubes):
    #
    # Compute distances to each cube
    #
    refdistances = np.sqrt(griddist[nx]*griddist[nx] + \
                           griddist[ny]*griddist[ny] + \
                           griddist[nz]*griddist[nz])
    for nspecies in range(nrspec):
        #
        # Check if sphere edge is inside the cube (can add "tolerance factors" here).
        #
        if  refdistances < inradius[nspecies] + 2.1*basecubesize \
        and refdistances > inradius[nspecies] - 2.1*basecubesize :
            refinearray1[nn] = 1
    #
    # First level disc refinements
    #
    if  griddist[nz] <  verticallimit0 \
    and griddist[nz] > -verticallimit0 \
    and refdistances <  radiallimit1   :
        refinearray1[nn] = 1
    #
    # Move coordinates
    #
    nx,ny,nz = movecoordinates(nx,ny,nz)
#
# Second refinements
#
refinearray2 = np.zeros(nbasecubes + 8*np.size(np.where(refinearray1 > 0)))
nx,ny,nz     = 0,0,0
counter      = 0
for nn in range(nbasecubes):
    #
    # Add 8 zeros for each 1 in refinement array.
    #
    if refinearray1[nn] == 1:
        refinearray2[nn+counter*8] = 1
        for nsmall in range(8):
            #
            # Check distance of refined cubes to centre (Sphere refinement)
            #
            refdistances = np.sqrt((griddist[nx]+gridcorrx1[nsmall])**2. + \
                                   (griddist[ny]+gridcorry1[nsmall])**2. + \
                                   (griddist[nz]+gridcorrz1[nsmall])**2.)
            for nspecies in range(nrspec):
                if  refdistances < inradius[nspecies] + 1.1*basecubesize \
                and refdistances > inradius[nspecies] - 1.1*basecubesize :
                    refinearray2[nn + counter*8 + nsmall+1] = 2
            #
            # Second level disc refinements
            #
            if  griddist[nz]+gridcorrz1[nsmall] <  verticallimit1 \
            and griddist[nz]+gridcorrz1[nsmall] > -verticallimit1 \
            and                    refdistances <    radiallimit2 :
                refinearray2[nn + counter*8 + nsmall+1] = 2
        counter += 1
    #
    # Move basegrid xyz-coordinates
    #
    nx,ny,nz = movecoordinates(nx,ny,nz)
#
# Third refinements
#
refinearray3      = np.zeros(nbasecubes + 8*np.size(np.where(refinearray2 > 0)))
nx,ny,nz          = 0,0,0
counter1,counter2 = 0,0
for nn in range(nbasecubes):
    if refinearray1[nn] == 1:
        refinearray3[nn+counter1*8+counter2*8] = 1
        for nsmall in range(8):
            if refinearray2[nn+counter1*8+nsmall+1] == 2:
                refinearray3[nn+counter1*8+counter2*8+nsmall+1] = 2
                for nsmall2 in range(8):
                    #
                    # Spherical refinements
                    #
                    refdistances = np.sqrt(\
                   (griddist[nx]+gridcorrx1[nsmall]+gridcorrx2[nsmall2])**2. + \
                   (griddist[ny]+gridcorry1[nsmall]+gridcorry2[nsmall2])**2. + \
                   (griddist[nz]+gridcorrz1[nsmall]+gridcorrz2[nsmall2])**2.)
                    for nspecies in range(nrspec):
                        if  refdistances < inradius[nspecies] + 0.6*basecubesize \
                        and refdistances > inradius[nspecies] - 0.6*basecubesize:
                            refinearray3[nn + counter1*8+counter2*8 + \
                                                           nsmall+nsmall2+2] = 3
                    #
                    # Third level disc refinements
                    #
                    if  griddist[nz] + gridcorrz1[nsmall] + gridcorrz2[nsmall2]\
                                                             <  verticallimit2 \
                    and griddist[nz] + gridcorrz1[nsmall] + gridcorrz2[nsmall2]\
                                                             > -verticallimit2 \
                    and                         refdistances <    radiallimit3 :
                        refinearray3[nn+counter1*8+counter2*8+nsmall+1+nsmall2+1] = 3
                counter2 += 1
        counter1 += 1
    #
    # Move basegrid xyz-coordinates
    #
    nx,ny,nz = movecoordinates(nx,ny,nz)
#
# Fourth and final refinements, only in the disc.
#
refinearray4               = np.zeros(nbasecubes + 8*np.size(np.where(refinearray3 > 0)))
nx,ny,nz                   = 0,0,0
counter1,counter2,counter3 = 0,0,0
for nn in range(nbasecubes):
    if refinearray1[nn] == 1:
        refinearray4[nn+counter1*8+counter2*8+counter3*8] = 1
        for nsmall in range(8):
            if refinearray2[nn + counter1*8 + nsmall+1] == 2:
                refinearray4[nn + 8*(counter1+counter2+counter3) + nsmall+1] = 2
                for nsmall2 in range(8):
                    if refinearray3[nn + 8*(counter1+counter2) + \
                                                         nsmall+nsmall2+2] == 3:
                        refinearray4[nn + 8*(counter1+counter2+counter3) + \
                                                           nsmall+nsmall2+2] = 3
                        for nsmall3 in range(8):
                            refdistances = np.sqrt(\
(griddist[nx]+gridcorrx1[nsmall]+gridcorrx2[nsmall2]+gridcorrx3[nsmall3])**2. + \
(griddist[ny]+gridcorry1[nsmall]+gridcorry2[nsmall2]+gridcorry3[nsmall3])**2. + \
(griddist[nz]+gridcorrz1[nsmall]+gridcorrz2[nsmall2]+gridcorrz3[nsmall3])**2.)
                            #
                            # Fourth level disc refinements
                            #
                            if  griddist[nz]        + gridcorrz1[nsmall ] + \
                                gridcorrz2[nsmall2] + gridcorrz3[nsmall3]   \
                                                          <  verticallimit3 \
                            and griddist[nz]        + gridcorrz1[nsmall ] + \
                                gridcorrz2[nsmall2] + gridcorrz3[nsmall3]   \
                                                          > -verticallimit3 \
                            and              refdistances <    radiallimit4:
                                refinearray4[\
nn + 8*(counter1+counter2+counter3) + nsmall+nsmall2+nsmall3+3] = 4
                        counter3 = counter3 + 1
                counter2 = counter2 + 1
        counter1 = counter1 + 1
    #
    # Move basegrid xyz-coordinates
    #
    nx,ny,nz = movecoordinates(nx,ny,nz)
#
# Define the final refined array for the third level.
#
refinearraysave = np.zeros(nbasecubes + 8*np.size(np.where(refinearray4 > 0)))
counter         = 0
for nn in range(np.size(refinearray4)):
    if refinearray4[nn] == 1:
        refinearraysave[nn + counter] = 1
    if refinearray4[nn] == 2:
        refinearraysave[nn + counter] = 1
    if refinearray4[nn] == 3:
        refinearraysave[nn + counter] = 1
    if refinearray4[nn] == 4:
        refinearraysave[nn + counter] = 1
        counter += 8
#
# Print amr_grid ============================================================= #
# --------------
# 1. The header
#
print("Writing amr_grid.inp")
nbranch = int(np.size(refinearraysave))
nleafs  = int(nbranch - (nbranch - nxyz*nxyz*nxyz) / 8)
f       = open('amr_grid.inp', 'w')
f.writelines(['1\n\n',\
              '1\n',\
              '0\n',\
              '0\n\n',\
              '1 1 1\n',\
              str(nxyz),' ',str(nxyz),' ',str(nxyz),'\n\n',\
              str(nrefines),' ',str(nleafs),' ',str(nbranch),'\n\n'])
#
# 2. The courners of the grid
#
for N in range(3):
    for nn in range(nxyz+1):
        f.writelines([str(gridcourners[nn]).strip('[]'),'\n'])
    f.writelines(['\n'])
f.writelines(['\n'])
#
# 3. The refinements
#
refinearraysave = refinearraysave.astype(int)
for nn in range(nbranch):
    f.writelines([str(refinearraysave[nn]).strip('[]'),'\n'])
f.close()
#
# Create density file ---------------------------------------------------------
# -------------------
#
# 1: Disc density file
#
# 1. Compute distances to each gridcell, including refined grid cells, daughter 
#    cells. Does it for disc and sphere parts simultaneously
#    Also adds disc densities, only for nspecies = 2 
#    (ie [1], since species nr 1 is [0])
#
print("Adding densities to the grid")
nspecies                               = 1
densitymatrix                          = np.zeros(nrspec*nleafs)
griddistances                          = np.zeros(nleafs)
nbig,nx,ny,nz                          = 0,0,0,0
counter1,counter2,counter3,disccounter = 0,0,0,0
for nn in range(np.size(refinearray1)):
    if refinearray1[nn] == 0:
        griddistances[nbig] = np.sqrt(griddist[nx]**2. + \
                                      griddist[ny]**2. + \
                                      griddist[nz]**2.)
        #
        # Add densities to base level (distances are to middle of cells why I
        # subtract half the cell vertically, done at every level)
        #
        if  griddist[nz] <  verticallimit0 - 0.5*basecubesize \
        and griddist[nz] > -verticallimit0 + 0.5*basecubesize \
        and griddistances[nbig] >            radiallimit1 \
        and griddistances[nbig] <     outradius[nspecies] :
            densitymatrix[nspecies*nleafs + nbig] = \
disczerodensity * (griddistances[nbig]/inradius[nspecies])**denspower[nspecies]
            disccounter += 1
        #
        # Move base coordinates
        #
        nx,ny,nz = movecoordinates(nx,ny,nz)
        #
        # Move mothergrid-index
        #
        nbig += 1
    #
    # Check for 1st level
    #
    if refinearray1[nn] == 1:
        #
        # Loop through 1st daughter cells (8 of them).
        #
        for nsmall in range(8):
            if refinearray4[nn + 8*(counter1+counter2+counter3) + nsmall+1] == 0:
                griddistances[nbig] = np.sqrt( \
                                     (griddist[nx] + gridcorrx1[nsmall])**2. + \
                                     (griddist[ny] + gridcorry1[nsmall])**2. + \
                                     (griddist[nz] + gridcorrz1[nsmall])**2.)
                #
                # Add densities to first level refinements
                #
                if  griddist[nz] <  verticallimit1 - 0.5*smallcubesize1 \
                and griddist[nz] > -verticallimit1 + 0.5*smallcubesize1 \
                and griddistances[nbig] >                  radiallimit2 \
                and griddistances[nbig] <                  radiallimit1 :
                    densitymatrix[nspecies*nleafs + nbig] = \
disczerodensity * (griddistances[nbig]/inradius[nspecies])**denspower[nspecies]
                    disccounter += 1
                #
                # Move mothergrid-index
                #
                nbig += 1
            #
            # Check for 2nd level
            #
            if refinearray4[nn + 8*(counter1+counter2+counter3) + nsmall+1] == 2:
                #
                # Loop through 2nd daughter cells (8 of them).
                #
                for nsmall2 in range(8):
                    if refinearray4[nn+8*(counter1+counter2+counter3)+nsmall+nsmall2+2]==0:
                        griddistances[nbig] = np.sqrt( \
               (griddist[nx] + gridcorrx1[nsmall] + gridcorrx2[nsmall2])**2. + \
               (griddist[ny] + gridcorry1[nsmall] + gridcorry2[nsmall2])**2. + \
               (griddist[nz] + gridcorrz1[nsmall] + gridcorrz2[nsmall2])**2.)
                        #
                        # Add densities to second level refinments
                        #
                        if  griddist[nz]+ gridcorrz1[nsmall] <  verticallimit2 - 0.5*smallcubesize2\
                        and griddist[nz]+ gridcorrz1[nsmall] > -verticallimit2 + 0.5*smallcubesize2\
                        and griddistances[nbig] > radiallimit3 \
                        and griddistances[nbig] < radiallimit2 :
                            densitymatrix[nspecies*nleafs + nbig] = \
disczerodensity * (griddistances[nbig]/inradius[nspecies])**denspower[nspecies]
                            disccounter += 1
                        #
                        # Move mothergrid-index
                        #
                        nbig += 1
                    #
                    # Check for 3rd level
                    #
                    if refinearray4[nn + 8*(counter1+counter2+counter3) +\
                                                         nsmall+nsmall2+2] == 3:
                        #
                        # Loop through 3rd daughter cells
                        #
                        for nsmall3 in range(8):
                            if refinearray4[nn + 8*(counter1+counter2+counter3)\
                                               + nsmall+nsmall2+nsmall3+3] == 0:
                                griddistances[nbig] = np.sqrt( \
(griddist[nx] + gridcorrx1[nsmall] + gridcorrx2[nsmall2] + gridcorrx3[nsmall3])**2. + \
(griddist[ny] + gridcorry1[nsmall] + gridcorry2[nsmall2] + gridcorry3[nsmall3])**2. + \
(griddist[nz] + gridcorrz1[nsmall] + gridcorrz2[nsmall2] + gridcorrz3[nsmall3])**2.)
                                #
                                # Add densities to third level refinments
                                #
                                if  griddist[nz]+ gridcorrz1[nsmall ] + \
                                    gridcorrz2[nsmall2] <  verticallimit3 - 0.5*smallcubesize3\
                                and griddist[nz]+ gridcorrz1[nsmall ] + \
                                    gridcorrz2[nsmall2] > -verticallimit3 + 0.5*smallcubesize3\
                                and griddistances[nbig] > radiallimit4 \
                                and griddistances[nbig] < radiallimit3 :
                                    densitymatrix[nspecies*nleafs + nbig] = \
disczerodensity * (griddistances[nbig]/inradius[nspecies])**denspower[nspecies]
                                    disccounter += 1
                                #
                                # Move mothergrid-index
                                #
                                nbig += 1
                            #
                            # Check for 4th level
                            #
                            if refinearray4[nn + 8*(counter1+counter2+counter3)\
                                               + nsmall+nsmall2+nsmall3+3] == 4:
                                for nsmall4 in range(8):
                                    #
                                    # Save radial distance of each cell
                                    #
                                    griddistances[nbig] = np.sqrt( \
              (griddist[nx] + gridcorrx1[nsmall ] + gridcorrx2[nsmall2] +      \
                              gridcorrx3[nsmall3] + gridcorrx4[nsmall4])**2. + \
              (griddist[ny] + gridcorry1[nsmall ] + gridcorry2[nsmall2] +      \
                              gridcorry3[nsmall3] + gridcorry4[nsmall4])**2. + \
              (griddist[nz] + gridcorrz1[nsmall ] + gridcorrz2[nsmall2] +      \
                              gridcorrz3[nsmall3] + gridcorrz4[nsmall4])**2.)
                                    #
                                    # Add disc's density to fourth refinements
                                    #
                                    if  griddist[nz]        + gridcorrz1[nsmall ]\
                                      + gridcorrz2[nsmall2] + gridcorrz3[nsmall3]\
                                           <  verticallimit4 - 0.5*smallcubesize4\
                                    and griddist[nz]        + gridcorrz1[nsmall ]\
                                      + gridcorrz2[nsmall2] + gridcorrz3[nsmall3]\
                                           > -verticallimit4 + 0.5*smallcubesize4\
                                    and griddistances[nbig] > inradius[nspecies] \
                                    and griddistances[nbig] <       radiallimit4 :
                                        densitymatrix[nspecies*nleafs + nbig] = \
disczerodensity * (griddistances[nbig]/inradius[nspecies])**denspower[nspecies]
                                        disccounter += 1
                                    #
                                    # Move mothergrid-index
                                    #
                                    nbig += 1
                        counter3 += 1
                counter2 += 1
        counter1 += 1
        #
        # Move base coordinates
        #
        nx,ny,nz = movecoordinates(nx,ny,nz)
print(f"Number of grid cells for the disc: {disccounter}")
#
# 2: Add in spherical densities in the densitymatrix, (species1 -> nspecies = 0)
#
# No need to recalculate the griddistances again.
# Compute densities as usual, according to equation defined in the beginning.
# and print file.
#
print('Writing dust_density.inp')
f = open('dust_density.inp', 'w')
f.writelines(['1\n',\
              str(nleafs),'\n',\
              str(nrspec),'\n'])
for nspecies in range(nrspec):
    for nn in range(nleafs):
        #
        # Set to zero outside the radial range.
        #
        if griddistances[nn] < inradius[nspecies] \
        or griddistances[nn] > outradius[nspecies]:
            densitymatrix[nspecies*nleafs + nn] = 0.
        #
        # Add the sphere densities to the densities
        #
        elif nspecies == 0:
            densitymatrix[nspecies*nleafs + nn] = \
                densitymatrix[nspecies*nleafs + nn] + spherezerodensity * \
               (griddistances[nn]/inradius[nspecies])**denspower[nspecies]
        #
        # Write the density file.
        #
        f.writelines([str(densitymatrix[nspecies*nleafs + nn]),'\n'])
f.close()
#
# Save array with grid distances in cm:
#
if savegrid == "y":
    print('Writing grid_distances.dat')
    np.savetxt('grid_distances.dat',griddistances,fmt='%f')
#
# Plot to check grid (Plotting takes time. Change nxyz > 1 to < 1 to save time!)
#
averageinradius  = np.mean(inradius)/AUcm
#
if plotgrid == "y":
    print("Plotting grid and densities")
    if nrspec == 2:
        plt.subplot(1,2,1)
        #
        # Disc densities
        #
        plt.plot(griddistances/AUcm,densitymatrix[nleafs:nrspec*nleafs],'ro',markersize=1, label='Disc')
        plt.plot([inradius[1]/AUcm,inradius[1]/AUcm],[0,np.max(densitymatrix)],'r')
        #
        # Sphere densities
        #
        plt.plot(griddistances/AUcm,densitymatrix[0:nleafs],'bo',markersize=1, label='Sphere')
        plt.plot([inradius[0]/AUcm,inradius[0]/AUcm],[0,np.max(densitymatrix)],'b:')
        #
        plt.legend()
        #
        # Density limits to work with.
        #
        plt.plot(averageinradius,0.1*totalzerodensity,'go',markersize=4)
        plt.plot(averageinradius,0.3*totalzerodensity,'go',markersize=4)
        plt.plot(averageinradius,0.5*totalzerodensity,'go',markersize=4)
        plt.plot(averageinradius,0.7*totalzerodensity,'go',markersize=4)
        plt.plot(averageinradius,0.9*totalzerodensity,'go',markersize=4)
        #
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'Radial distance to centrum of grid (AU)')
        plt.ylabel(r'Dust density (g/cm$^3$)')
        #
        # And to check disc thickness
        #
        plt.subplot(1,2,2)
        plt.plot([ inradius[1]/AUcm,radiallimit4/AUcm] , [verticallimit4/AUcm,verticallimit4/AUcm],'k')
        plt.plot([radiallimit4/AUcm,radiallimit3/AUcm] , [verticallimit3/AUcm,verticallimit3/AUcm],'k')
        plt.plot([radiallimit3/AUcm,radiallimit2/AUcm] , [verticallimit2/AUcm,verticallimit2/AUcm],'k')
        plt.plot([radiallimit2/AUcm,radiallimit1/AUcm] , [verticallimit1/AUcm,verticallimit1/AUcm],'k')
        plt.plot([radiallimit1/AUcm,outradius[0]/AUcm] , [verticallimit0/AUcm,verticallimit0/AUcm],'k')
        plt.plot([0,outradius[1]/AUcm],[0,0.5*discflare*outradius[1]/AUcm],'k:')
        plt.xlabel(r'Radial distance to centrum of grid (AU)')
        plt.ylabel(r'Disc height, $h(r)$ (from mid-disc-plane, AU)')
#
# Write info file for reference.
#
print("Writing grid_info.txt")
f = open('grid_info.txt', 'w')
f.writelines([\
 '               Sphere between: ',str(inradius[0]/AUcm),' AU to ',str(outradius[0]/AUcm),' AU\n'\
,'                 Disc between: ',str(inradius[1]/AUcm),' AU to ',str(outradius[1]/AUcm),' AU\n'\
,'Approx. disc flaring relation: ',str(discflare),' * r\n'\
,'            Sphere mass ratio: ',str(sphereratio),'\n'\
,'                Input mass is: ',str(totaldustmass/1.989e33/spheremasscorr),' Msol\n'\
,'Sphere, species1 density distr.: propto r^',str(denspower[0]),'\n'\
,'  Disc, species2 density distr.: propto r^',str(denspower[1]),'\n'\
,'                 N disc cells: ',str(disccounter),'\n'\
,'\n'\
,'Density percentage\n'\
,'sphere / (disc + sphere) = ',str(spherezerodensity / (disczerodensity + spherezerodensity)),'\n'\
,'disc   / (disc + sphere) = ',str(disczerodensity / (disczerodensity + spherezerodensity)),'\n'\
,'\n'\
,'     Final smallest grid size: ',str(smallcubesize4/AUcm),'AU\n'\
,'      Third refined grid size: ',str(smallcubesize3/AUcm),'AU\n'\
,'     Second refined grid size: ',str(smallcubesize2/AUcm),'AU\n'\
,'      First refined grid size: ',str(smallcubesize1/AUcm),'AU\n'\
,'            Largest grid size: ',str(basecubesize/AUcm),'AU\n'\
,'Disc half thicknesses:\n'\
,' Inside Radius =',str(radiallimit4/AUcm),'AU, Thickness:',str(verticallimit4/AUcm),'AU, Number of 4th refined cells:',str(numbcells4),'\n'\
,' Inside Radius =',str(radiallimit3/AUcm),'AU, Thickness:',str(verticallimit3/AUcm),'AU, Number of 3rd refined cells:',str(numbcells3),'\n'\
,' Inside Radius =',str(radiallimit2/AUcm),'AU, Thickness:',str(verticallimit2/AUcm),'AU, Number of 2nd refined cells:',str(numbcells2),'\n'\
,' Inside Radius =',str(radiallimit1/AUcm),'AU, Thickness:',str(verticallimit1/AUcm),'AU, Number of 1st refined cells:',str(numbcells1),'\n'\
,'Outside Radius =',str(radiallimit1/AUcm),'AU, Thickness:',str(verticallimit0/AUcm),'AU, Number of base cells       ',str(numbcells0)])
f.close()
print('All done!')
os.system('spd-say "Moo"')
plt.show()
