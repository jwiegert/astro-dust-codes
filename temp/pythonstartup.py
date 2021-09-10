# Useful libraries ============================================================#
# ----------------                                                             #
import os
# Numpy: ----------------------------------------------------------------------#
import numpy as np
from numpy.matlib import rand,zeros,ones,empty,eye
# Scipy:-----------------------------------------------------------------------#
import scipy as s
from scipy.integrate import quad
from scipy import ndimage
from scipy.interpolate import interp1d
#from scipy import optimize
## Astronomy packages ----------------------------------------------------------#
#import pyfits as pf
#from astropy import *
#from astropy import units as u
#import radmc3dPy as r3d
#import astropysics as ast
## Some extra commands ---------------------------------------------------------#
import math
import random
# Plotting tool ---------------------------------------------------------------#
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['serif']})
rc('text', usetex=True)
# Useful Constants ============================================================#
# ----------------                                                             #
from constants import *
#==============================================================================#
