# kinetics-py
# Phase transformation kinetics

### Introduction

This folder holds the python code to fit kinetics data and generate the figures of the following publication
> C. Langrand, D. Andrault, S. Durand, Z. Konôpková, N. Hilairet, C. Thomas & S. Merkel (2019)
> Kinetics and detectability of the bridgmanite to post-perovskite transformation in the Earth's D′′ layer 
> Nat. Commun. , 10 [doi: 10.1038/s41467-019-13482-x] 
> [https://www.nature.com/articles/s41467-019-13482-x](https://www.nature.com/articles/s41467-019-13482-x)

The python program will read kinetics data from a file with a list of temperatures, pressures, and transformation timescales, phase proportions, and change of unit cell volumes and fit kinetics law based on
* a nucleation and growth model
* a shear transformation model

The fitted parameters will be written to the screen. The program will also generate figures such as those of Fig. 3 in the publication of Langrand et al.

### List of files
* hernlundcold.dat: the cold geotherm of Hernlund et al, 2005
* hernlundhot.dat: the hot geotherm of Hernlund et al, 2005
* hernlundwarm.dat: the mean geotherm of Hernlund et al, 2005
* kinetics-data.dat: the kinetics data of Langrand et al, 2019
* no-kinetics-data.dat: the P/T points for which Langrand et al, 2019, did not observe any transformation
* kinetics-fits.py: the python code

### Example

To fit kinetics and generate plots assuming a Pv/pPv phase boundary given by a reference P/T point at 128 GPa and 3300 K, with a Clapeyron slope of 8.5 MPa/K, type
> /kinetics-fits.py -P 128 -T 3300 -c 8.5e-3 -n no-kinetics-data.dat kinetics-data.dat

This will output the fit results to the screen, plot the kinetics models, and generate files for the corresponding figures, in SVG and PNG format.

Generate figures will look like the one below
![Example results for a nucleation and growth model](https://github.com/smerkel/kinetics-py/blob/master/example-nuclGrowthModel.png)

### Full instructions

This is the ouput of `./kinetics-fits.py -h`:

    usage: kinetics-fits.py [options] datafile
    
    Fit of a shear or nucleation and growth model for the perovskite to post-
    perovskite transformation. Adapted from C. Langrand PhD thesis.
    
    positional arguments:
      datafile              Data file (required)
    
    optional arguments:
      -h, --help            show this help message and exit
      -P PE, --Pe PE        Pressure of point on pv / (pv+ppv) phase boundary in
                            GPa. Default is 128.0
      -T TE, --Te TE        Temperature of point on pv / (pv+ppv) phase in K.
                            Default is 3300.0
      -c CLAPEYRON, --clapeyron CLAPEYRON
                            Clapeyron's slope in GPa/K. Default is 0.0067
      -f FACTOR, --factor FACTOR
                            Multiply times by 10^x in plotting. Default is 0.0
      -n NONEOBSERVATION, --noneobservation NONEOBSERVATION
                            File with negative observations. Default is None

