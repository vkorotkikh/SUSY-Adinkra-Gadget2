# Adinkra Gadget II Calculation Script
Description: This repository contains the Python code for calculating all the
"Gadget II" values over all the possible BC4 space Adinkras. Given 36,864 Adinkras
which are unique sets of four n=4 matrices, there are 36,864^2 possible Gadget II
values. 

## DEV STATUS
```
Use Python 3.5 or newer for testing or calculating Gadget II values. Found an
unknown bug using Python 2.7 that produces incorrect Gadget II values.
Updating pieces of code for public release.
```


## Important details
Developed and tested using Python 3 (using Anaconda Python 3.5 as of README update).
```
$ Requires numpy library
$ Developer/creator tested on Mac OSX.
```

## Getting Started
To execute this tool, run the following command in Terminal
```
$ python -u run_adinkra_calc.py
```
or if you have both python 2 and 3, specify python 3 version.
```
$ python3 run_adinkra_calc.py
```


## General Overview
### Primary script
adinkra_tetrad_calc.py - Calculates the 36,864 ordered BC4-based
adinkras with 4 colours, 4 open and 4 closed nodes. In this case each Adinkra is
a tetrad of 4 L matrices.
### Code for Gadget calculation
fx_mpgadgets.py - Contains the multiprocessing enabling code for fast calculating
and saving all of the Gadget II results.

## Author

-- Vadim Korotkikh --
