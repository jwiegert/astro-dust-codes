# Readme: astro-dust-codes

This repository contains a collection of Python scripts used for my research on circumstellar dust. These scripts are primarily to create input data for radiative transfer simulations with [RADMC-3D](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/) currently at v2.0 (a MC-based raytracing based radiative transfer simulation package written by C.P. Dullemond at Heidelberg), and for visualising input and analysing output spectra and images. Most often I use R3D for simulating SEDs at a wide wavelength range.

Between 2011 and 2016, I mostly worked with dust of cirumstellar debris discs around sun-like stars, for these I used a restricted N-body simulation I wrote in C to create my dust models. Since 2018, my work has been on dust formed in the envelopes of AGB-star.

For a list of publications, please see my [portfolio-readme](https://github.com/jwiegert/portfolio-joachim-wiegert/blob/main/README.md).

To Be Added:

- Scripts and description.

- Instructions on how to use the different files, how to create 3D visulisations of the models.

- Required packages. (Numpy, SciPy, Matplotlib).

## List of scripts

- [constants.py](https://github.com/jwiegert/astro-dust-codes/blob/main/constants.py)

A list of useful constants in astronomy etc.

- [inputdata.py](https://github.com/jwiegert/astro-dust-codes/blob/main/inputdata.py)

Settings and properties used in these scripts to create input data for RADMC-3D. Lists stellar propertes, dust grain properties (used for RADMC-3D's included Mie-code to create opacity files), grid properties, dust distribution properties. (It's a bit messy now but it works.)

- [createsphereflaring4stepsdisccombineNspecies4refined.py](https://github.com/jwiegert/astro-dust-codes/blob/main/createsphereflaring4stepsdisccombineNspecies4refined.py)

File name describes function of the script. This creates a grid and density file for RADMC-3D where the dust is distributed in a sphere and a disc (with different inner and outer radii), powerlaw density distributions (density proportional to radius^-k). The disc is approximating a flared disc. The grid is octree refined to four levels with finer details in the centre of the grid. It includes up to two dust species separated into the sphere and disc distributions. Can also be evenly distributed. All properties and settings are changed in the [inputdata.py](https://github.com/jwiegert/astro-dust-codes/blob/main/inputdata.py) file.

- [plotimageout.pu]((https://github.com/jwiegert/astro-dust-codes/blob/main/plotimageout.py)

Plots image output from RADMC-3D. Filenames need to be changed manually in script and flux limits may be changed manually also.

## Additional instructions

- How to visualise your input data.

RADMC-3D can create *.vtk (Visualization Toolkit format file) files using the command

> radmc3d vtk_dust_density 1

for each of the densities. There are a number of programs to open and visualise these. The open source program [paraview](https://www.paraview.org/) is been useful for this.

