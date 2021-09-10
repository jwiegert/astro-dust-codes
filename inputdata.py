# Input data used to create *.inp-files necessary to run radmc-3d.             #
#                                                                              #
# Required input files to run radmc-3d are:                                    #
#   amr_grid.inp                                                               #
#   dust_density.inp                                                           #
#   dustkappa_*.inp   FOR each specie                                          #
#   dustopac.inp                                                               #
#   radmc3d.inp                                                                #
#   stars.inp                                                                  #
#   wavelength_micron.inp                                                      #
# Then you run radmc3d mctherm to create                                       #
#   dust_temperature.dat                                                       #
#                                                                              #
#==============================================================================#
#
AUcm = 14959800000000.0
#
# Stellar properties. ======================================================== #
#                                                                              #
distance         = 8200.0             # In pc
distanceerror    =  820.0             # error in pc
nstar            = 1                  # Number of stars.
rstar            = 323.*6.955e10      # Radius of each star in cm.
mstar            = 2.*1.989e33        # Mass of each star in gram.
coordinatesstars = [0.,0.,0.]         # [x,y,z] coordinate of each star in cm.
                                      # (Add rows for more than one star.)
startemperature  = 2600.              # Black body temperature of each star.
luminosity       = 4110. * 3.846e26   # Solar luminosities
#
# Dust distribution. ========================================================= #
#                                                                              #
nrspec            = 2                         # Nr of species
#
# Dust temperature at the inner radius set to 1000K
# about 13.8 Rstar.
# Router is about 3700 Rinner (or about 20K)
#
Tin               = 1000.
Tout              =   20.
#
# ---------------------------------------------------------------------------- #
#
# Sphere Mass corrections (20 to 4000 AU distribution)
# 13.8Rstar ~ 20.7 AU ~ 20AU.
#
spheremasscorr    =  1.000992076313
inradius          = [      20.*AUcm,     20.*AUcm]  # Inner radius.
#
#outradius         = [    4000.*AUcm ,  4000.*AUcm]  # Outer radius.
#discmasscorr      = 14.490763389451 # 4000au, 10deg
#sphereratio       =  0.3
#discmasscorr      =   # 4000au, 30deg
#discmasscorr      = 2.2498809358876 # 4000au, 45deg
#
#outradius         = [    4000.*AUcm ,  3200.*AUcm]  # Outer radius.
#discmasscorr      = 14.600828245876 # 3200au, 10deg
#sphereratio       =  0.1
#discmasscorr      =   # 3200au, 30deg
#discmasscorr      =  2.206789304898 # 3200au, 45deg
#
#outradius         = [    4000.*AUcm ,  2400.*AUcm]  # Outer radius.
#discmasscorr      = 16.582646514113 # 2400au, 10deg
#discmasscorr      =   # 2400au, 30deg
#discmasscorr      =  2.341372424164 # 2400au, 45deg
#
#sphereratio       =  0.3 # Fraction of sphere mass - gives approx 3% of density in sphere.
#
# DOWN FROM HERE
#
outradius         = [    4000.*AUcm ,  1600.*AUcm]  # Outer radius.
discmasscorr      =  6.331768198315633 # 1600au, 10deg
sphereratio       = 0.993 # Disc density: 10%
#sphereratio       = 0.9738 # Disc density: 30%
#sphereratio       = 0.941 # Disc density: 50%
#sphereratio       = 0.8725 # Disc density: 70%
#sphereratio       = 0.638 # Disc density: 90%
#
#discmasscorr      =  1.312249718559036 # 1600au, 45deg
#sphereratio       = 0.968 # Disc density: 10%
#sphereratio       = 0.885 # Disc density: 30%
#sphereratio       = 0.768 # Disc density: 50%
#sphereratio       = 0.587 # Disc density: 70%
#sphereratio       =  0.27 # Disc density: 90%
#
#
#outradius         = [    4000.*AUcm ,   800.*AUcm]  # Outer radius.
#discmasscorr      = 5.64628712181288 #  800au, 10deg
#sphereratio       = 0.996 # Disc density: 10% 
#sphereratio       = 0.9854 # Disc density: 30%
#sphereratio       = 0.9665 # Disc density: 50%
#sphereratio       = 0.925 # Disc density: 70%
#sphereratio       =  0.76 # Disc density: 90% 
#
#discmasscorr      =  1.27819628121520 #  800au, 45deg
#sphereratio       =  0.9834 # Disc density: 10%
#sphereratio       = 0.9385 # Disc density: 30%
#sphereratio       = 0.867 # Disc density: 50%
#sphereratio       = 0.737 # Disc density: 70%
#sphereratio       =  0.42 # Disc density: 90%
#
#
#outradius         = [    4000.*AUcm ,   200.*AUcm]  # Outer radius.
#discmasscorr      =  3.29668678654542 #  200au, 10deg
#sphereratio       = 0.9985 # Disc density: 10%
#sphereratio       = 0.99415 # Disc density: 30%
#sphereratio       = 0.9865 # Disc density: 50%
#sphereratio       = 0.969 # Disc density: 70%
#sphereratio       = 0.89 # Disc density: 90%
#
#discmasscorr      =  1.09169140023270 #  200au, 45deg
#sphereratio       =  0.9955 # Disc density: 10%
#sphereratio       = 0.9825 # Disc density: 30%
#sphereratio       = 0.9602 # Disc density: 50%
#sphereratio       = 0.912 # Disc density: 70%
#sphereratio       =  0.727 # Disc density: 90%
#
#
#outradius         = [    4000.*AUcm ,   100.*AUcm]  # Outer radius.
#discmasscorr      =  2.10652021246658 #  100au, 10deg
#sphereratio       =  0.999 # Disc density: 10%
#sphereratio       = 0.99592 # Disc density: 30%
#sphereratio       = 0.99055 # Disc density: 50%
#sphereratio       = 0.9783 # Disc density: 70%
#sphereratio       =  0.921 # Disc density: 90%
#
#discmasscorr      =  1.00204596991819 #  100au, 45deg
#sphereratio       = 0.9979 # Disc density: 10%
#sphereratio       = 0.9915 # Disc density: 30%
#sphereratio       = 0.9803 # Disc density: 50%
#sphereratio       = 0.9553 # Disc density: 70%
#sphereratio       =  0.847 # Disc density: 90%
#
denspower         = [-2.,-2.]                         # Density power law.
totaldustmass     = 1.e-3 * 1.989e33 * spheremasscorr # Total mass in gram
#discflare         = 2.* 0.087489   #  5 deg
discflare         = 2.*0.17633     # 10deg; For discs, h(r) = discflare * r
#discflare         = 2.*1.0         # 45deg; For discs, h(r) = discflare * r
#
# tan(flareangle) = 2.*discflare
#  5 deg = 0.087489
# 10 deg = 0.17633
# 30 deg = 0.57735
# 45 deg = 1.0
#
#discwidth         = 0.5 * inradius[0]                  # For DISCS, half width of disc (half h(r)). - olddiscs
#
# Notes about
# * createsphereflaringdisccombineNspecies4refined.py
#
# These file puts the Spherical distribution as species number1.
# The disc distribution is put as species number 2.
#
# Wavelength-grid.============================================================ #
#
lambdazero       = 0.1   # Wavelength range in micrometers
lambdafinal      = 4000. #
nlambda          =  500  # Number of wavelength points
logarithmicrange = 1     # Chose, 0 for linear range, 1 for logarithmic range
#
# Grain properties.=========================================================== #
#
# To use if you need to create new mass opacity files (kappa_abs, kappa_scat)
# with Mie-theory - using files included in the radmc3d-package.
#                                                                              #
#agrainmax        = 1e-4*100.   # Maximum size in cm
#na               = 9           # Number of grain size samples.
#agrainmin        = 2e-7*100.   # Minimum grain size in cm
#agrainpw         = -3.5        # Grain size power law: N propto agrain^agrainpw
#
# Optical constant file names (drop .lnk)
#
#optconst         = 'mg2sio4'    #  Mg(2) SiO(4)  [? g/ccm]
#optconst         = 'fe2sio4' #  Mg Fe SiO(4) [3.71 g/ccm]
#optcont          = 'al2o3porous' #
#optconst         = 'mg2sio4toy'    #  Mg(2) SiO(4)  [? g/ccm]
#matdens          = 2.65    # The material density in gram / cm^3
