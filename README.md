# SPyCFit

## Overview

**SPyCFIT** is the **S**imple **Py**thon **C**osmology **Fit**ter I am going to try and write to fit cosmological models to supernova data from various lightcurve fitters.

## Plan

Once finished the code will allow you to

* Read in low-z data from any known fitter (SALT2, SiFTO or SNOOPY)
* add high-z data
* fix any cosmology parameters or lightcurve nuisance parameters
* do a PyMinuit-based chi-squ minimisation over any free parameters by tweeking the instrinsic dispersion of the SN until $$$\chi^2_r = 1$$$.

## Progress

3/10/13 - Use argparse to get optional arguments

3/10/13 - amalgamate_salt2.py organises Salt2 data so I can start writing and testing the main code

10/10/13 - amalgamate_salt2.py is now a stand-alone programme.  This makes the main code much easier to write

10/10/13 - use setattr() to read the output of amalgamate SALT2 into the classes and create the snvals dict with entry per object.

10/10/13 - fixed some argparse weirdness.  I didn't understand why it was doing what it was, but it works fine now.

15/05/14 - this is fully-working for SALT2.4 data with LMFIT as the minimiser for the $$$\chi^2$$$ of the fit.  Full options are described below

## Maths

The code is designed to minimise the formula

$$\sum \frac{5\log \mathcal{D}_L(z,\theta) - (m-\mathcal{M})}{\sigma_s^2 + \sigma_i^{2} + \sigma_z^2 }$$

where $$$m-\mathcal{M}$$$ is the observed distance modulus where $$$m$$$ is the corrected B-band magnitude at maximum and $$$\mathcal{M} = M_B -5\log H_0 + 25$$$.  The errors are the statistical, intrinsic dispersion and the error on the redshifts.

For the luminosity distance, I used a Hubble-constant-less version of the forumula (eg the one at the bottom of the second page of Perlmutter et al. 1997, ApJ, 483, 565) and multiplied through by $$$H_0$$$ so that

$$\mathcal{D}_L = H_0 d_L$$

The free parameters in the LMFIT are $$$\alpha,\beta$$$ (from the light curve standardisation) and $$$\mathcal{M}$$$.  You adjust the value of $$$\sigma_i$$$ by hand until $$$\chi^2_r = 1$$$.

## Input File Format

The input file to spycfit.py should be from __1__ LC fitter only which is set by the style parameter.  Line 1 of the file must be the names of the columns.  The first column __must__ be the name of the supernova.  The rest of the columns can follow in any order, but for a given LC fitter they must contain:

**SALT2**

* Object with column 1 heading name
* Redshift and redshift error with respective column headings z and z_err
* B-band max date and error headed as bmax and bmax_err
* B-band maximum observed magnitude
* X0 and error as x0 and x0_err
* X1 and error as x1 and x1_err
* Colour and error named as colour and colour_err (note the u in colour!)

**SiFTO**

* Object with column 1 heading name
* edshift and redshift error with respective column headings z and z_err
* B-band max date and error headed as bmax and bmax_err
* B-band maximum observed magnitude
* Stretch and error as stretch and stretch_err
* Colour and error named as colour and colour_err (note the u in colour!)

**SNOOPY**

To do

## Options

The options are parsed using argsparse.

usage: spycfit.py [-h] [-s SIGMA] [-a ALPHA] [-b BETA] [-m OMEGA_M]
                  [-l OMEGA_L] [-f] [-z] [-e EXPANSION_RATE] [-c CUTS] [-w]
                  filename {salt2,sifto,snoopy}

Simply Python Cosmology Fitter

positional arguments:
  filename              Input file of SN data. Name in first column
  {salt2,sifto,snoopy}  Style of input data. NB Snoopy not implemented yet!

optional arguments:

  -h, --help            show this help message and exit
  
  -s SIGMA, --sigma SIGMA Instrinsic dispersion
  
  -a ALPHA, --alpha ALPHA Fix coeff for LC width
  
  -b BETA, --beta BETA  Fix coeff for LC colour
  
  -m OMEGA_M, --omega_m OMEGA_M Fix value of Omega_matter
  
  -l OMEGA_L, --omega_l OMEGA_L Fix value of Omega_lambda
  
  -f, --flat            Force a flat universe model
  
  -z, --highz           Add high-z data (not implemented)
  
  -e EXPANSION_RATE, --expansion_rate EXPANSION_RATE Value of the Hubble constant
  
  -c CUTS, --cuts CUTS  Apply quality cuts. No cuts is default
  
  -w, --writetofile     Write out plot to eps file
  
Syntax examples:

run spycfit.py salt24_050514.dat salt2  -f -m 0.3 -l 0.7 -s 0.14 -c betoule_salt2.dat -w

This is the best way to run the code.  It uses SALT24 data in a fixed, flat cosmology with the instrinsic dispersion also fixed.  The alpha and beta parameters are free.  Cuts from the Bertoule paper are used.  The figure is saved to eps.

run spycfit.py salt24_050514.dat salt2  -f -m 0.25 -l 0.75 -s 0.14 -a 0.1

A fixed flat cosmology and the alpha parameter fixed.

## Required Packages

* numpy
* matplotlib
* sys
* os
* argparse
* \_\_future__
* lmfit 
* scipy
* import datetime

