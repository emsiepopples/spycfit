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

class Supernova(object):
	
	
	name=''
	z = float('nan')
	z_err = float('nan')
	maxdate = float('nan')
	maxdate_err = float('nan')
	bmax = float('nan')
	bmax_err = float('nan')
	fitter = ''
	filepath = ''
	
	__doc__ = "Supernova Class"
	

class Salt2SN(Supernova):
	
	#initialise object as empty
	
	fitter = 'salt2'
	colour = float('nan')
	colour_err = float('nan')
	x0 = float('nan')
	x0_err = float('nan')
	x1 = float('nan')
	x1_err = float('nan') 
	
	__doc__ = "SALT2 subclass of Supernova Class"
	
	
	

		


class SiftoSN(Supernova):

	#initialise object as empty
	
	fitter = 'sifto'
	colour = float('nan')
	colour_err = float('nan')
	stretch = float('nan')
	stretch_err = float('nan')

	#will need a getdata method as for the SALT data
	
	__doc__ = "SiFTO subclass of Supernova Class"
	


def main(argv=None):

	parser = argparse.ArgumentParser(description='Simply Python Cosmology Fitter')
	parser.add_argument("filename", type=str, help="Input file of SN data. Name in first column")
	parser.add_argument("style", type=str, choices=['salt2','sifto', 'snoopy' ],help="Style of input data")
	parser.add_argument("sigma", type=float, help='Instrinsic dispersion')
	parser.add_argument("-a", "--alpha", type=float, help='Fix coeff for LC width')
	parser.add_argument("-b", "--beta", type=float, help='Fix coeff for LC colour')
	parser.add_argument("-m", "--omega_m", type=float, help='Fix value of Omega_matter')
	parser.add_argument("-l", "--omega_l", type=float, help='Fix value of Omega_lambda')
	parser.add_argument("-z", "--highz", action="store_true", help='Add high-z data')
	parser.add_argument("-e", "--expansion_rate", type=float, help='Value of the Hubble constant')
	
	args = vars(parser.parse_args())

#read in the datafile, getting column names from row 0 and save into a structured array

	f = open(args['filename'], 'r')
	col_names = f.readline().split()
	lines = f.readlines()
	lines = [x for x in lines if not x.startswith('#')]  #removes comments
	f.close()
	snvals = {}
	
	
	for ll in lines:
		
		dat = ll.split()
		
		if args['style'] == 'salt2':
			snvals[dat[0]] = Salt2SN()
		elif args['style'] == 'sifto':
			snvals[dat[0]] = SiftoSN()
		elif args['style'] == 'snoopy':
			snvals[dat[0]] = SnoopySN()
			
		
			
		for idx,col in enumerate(col_names): 
			
			if idx == 0: 
				setattr(snvals[dat[0]], col, dat[idx])
			else:
				setattr(snvals[dat[0]], col, float(dat[idx]))
		
		
		
#can now cycle of snvals.keys() for each object!!
#need to be able to use class or structured array on the high-z data as well.
#construct function and minuit call for minimisation.

#	for kk in snvals.keys():
#		print snvals[kk].colour

			
"""
def chi2(some arguments):
	chi2 = sum((data-model)**2./(sigma**2 + int_disp**2))
	
def import_highz():
	work out some way to import the high-z data to add to the minimisation	
"""

if __name__ == "__main__":
	main()