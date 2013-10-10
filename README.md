# SPyCoFit

## Overview

**SPyCoFIT** is the **S**imple **Py**thon **Co**smology **Fit**ter I am going to try and write to fit cosmological models to supernova data from various lightcurve fitters.

## Plan

Once finished the code will allow you to

* Read in low-z data from any known fitter (SALT2, SiFTO or SNOOPY)
* add high-z data
* fix any cosmology parameters or lightcurve nuisance parameters
* do a PyMinuit-based chi-squ minimisation over any free parameters by tweeking the instrinsic dispersion of the SN until the reduced chi-squ = 1

## Progress

3/10/13 - Use argparse to get optional arguments

3/10/13 - amalgamate_salt2.py organises Salt2 data so I can start writing and testing the main code

10/10/13 - amalgamate_salt2.py is now a stand-alone programme.  This makes the main code much easier to write

10/10/13 - use setattr() to read the output of amalgamate SALT2 into the classes and create the snvals dict with entry per object.

10/10/13 - fixed some argparse weirdness.  I didn't understand why it was doing what it was, but it works fine now.