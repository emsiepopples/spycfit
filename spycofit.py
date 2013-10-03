#/usr/bin/python

"""
SPyCoFit
Simple Python Cosmology Fitter written by Emma Walker

"""


import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import argparse
import __future__

import pdb

def main(argv=None):

	parser = argparse.ArgumentParser()
	parser.add_argument("filename", type=str, help="Input file of SN data")
	parser.add_argument("style", type=str, choices=['salt2','sifto', 'snoopy' ],help="Style of input data")
	parser.add_argument("sigma", type=float, help='Instrinsic dispersion')
	parser.add_argument("-a", "--alpha", type=float, help='Fix coeff for LC width')
	parser.add_argument("-b", "--beta", type=float, help='Fix coeff for LC colour')
	parser.add_argument("-m", "--omega_m", type=float, help='Fix value of Omega_matter')
	parser.add_argument("-l", "--omega_l", type=float, help='Fix value of Omega_lambda')
	parser.add_argument("-z", "--highz", action="store_true", help='Add high-z data')
	parser.add_argument("-h", "--hubble", action=float, help='Value of the Hubble constant')
	
	args = parser.parse_args()
#	print args()
	if args.omega_m:
		print("Fixed Omega_m: {0}".format(args.omega_m))
		
#next step: work out the best way to read the data in. Create class?
#will need a big array of data.
#need to be able to use class or structured array on the high-z data as well.
#construct function and minuit call for minimisation.
			
"""
def chi2(some arguments):
	chi2 = sum((data-model)**2./(sigma**2 + int_disp**2))
	
def import_highz():
	work out some way to import the high-z data to add to the minimisation	
"""

if __name__ == "__main__":
	sys.exit(main())