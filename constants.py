# Standard constants
#
# (i) Constants
#
c          = 2.998e8
kb         = 1.3806503e-23
hplanck    = 6.626068e-34
G          = 6.67384e-11
Gcgs       = 6.67384e-8
sigmaSB    = 5.670373e-8
sigmaSBcgs = 5.670373e-5
#
# (ii) Units
#
AU         = 1.49598e11
AUcm       = 1.49598e13
pc         = 30.857e15
eV         = 1.6021766e-19
asecsqster = 2.35044305391e-11 # 1 asec squared equals this in sterradian
#
# (iii) Properties
#
Mmoon      = 7.3477e22
Mearth     = 5.9736e24
Msol       = 1.989e30
Lsol       = 3.846e26
Rsol       = 6.955e8
Tsol       = 5778.0
loggsun    = 4.440
#
# (iv) Projects
#      alpha Centauri
MacenA     = 1.105 * Msol
MacenB     = 0.934 * Msol
eccacen    = 0.5179
acenmajor  = 23.684
#
#      Complex fields
#      - Hip 4148
Lh4148     = 0.297 * Lsol
Mh4148     = 0.74 * Msol
Rh4148     = 0.79 * Rsol
Teffh4148  = 4940.
Disth4148  = 12.22 * pc
Disth4148Error = 1.58 * pc
#
#      - Hip 13402
#
Lh13402 = 0.392 * Lsol
Mh13402 = 0.89 * Msol
#Rh13402 = 1.1 * Rsol
Rh13402 = 0.8 * Rsol
Teffh13402 = 5217.

Disth13402 = 10.35 * pc
Disth13402Error = 0.04 * pc

#      - Hip 14954
Lh14954A = 3.848 * Lsol
Mh14954A = 1.34 * Msol
Rh14954A = 1.82 * Rsol
Teffh14954A = 6187.

#Lh14954B = 0.0927 * Lsol  # From Stefan-Boltzmann's law
#Lh14954B = 0.009 * Lsol   # From flux density ratios and a total 0.011
Lh14954B = 0.05 * Lsol     # From tables on M dwarfs
Mh14954B = 0.554 * Msol    # From their orbits
#Rh14954B = 0.52 * Rsol    # From Fracassini et al 1988 (why did I have 0.55 before?)
Rh14954B = 0.55 * Rsol     # From stefan-boltzmann's law
#Teffh14954B = 4300.       # From Roell
Teffh14954B = 3700.        # From Tables

#Lh14954C = 0.0269 * Lsol  # From flux ratio C/B = 0.29 (Roell 201X)
#Lh14954C = 0.002 * Lsol   # From flux ratio C/B = 0.29 (Roell 201X) and a total of 0.011
Lh14954C = 0.01 * Lsol     # From tables on M dwarfs
Mh14954C = 0.336 * Msol    # From their orbits
#Rh14954C = 0.4227 * Rsol  # From SB-law and flux ratio.
Rh14954C = 0.31 * Rsol     # From SB-law
#Teffh14954C = 3600.       # From Roell
Teffh14954C = 3300.        # From tables

Disth14954 = 22.6 * pc
Disth14954Error = 0.1 * pc
