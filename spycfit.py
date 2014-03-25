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
from lmfit import minimize, Parameters, report_fit
import pdb
from scipy import integrate
from IPython import embed

class Supernova(object):
	
	def __init__(self, arr):
	
		self.name=arr[0]
		self.z = float(arr[1])
		self.z_err = float(arr[2])
		self.maxdate = float(arr[3])
		self.maxdate_err = float(arr[4])
		self.bmax = float(arr[5])
		self.bmax_err = float(arr[6])
		
	def integrand(self, omega_m, omega_l):
	
		return  1.0/(np.sqrt((1.+self.z)**2 * (1+ omega_m) - z*(2.+self.z)*omega_l))
		
	
	__doc__ = "Supernova Class"
	

class Salt2SN(Supernova):
	
	#initialise object as empty
	
	def __init__(self, arr):
		
		super(Supernova, self).__init__(arr)
		
		self.fitter = 'salt2'
		self.colour = float(arr[11])
		self.colour_err = float(arr[12])
		self.x0 = float(arr[7])
		self.x0_err = float(arr[8])
		self.x1 = float(arr[9])
		self.x1_err = float(arr[10]) 

	
	__doc__ = "SALT2 subclass of Supernova Class"
	


class SiftoSN(Supernova):
	
	def __init__(self,arr):
		
		super(Supernova, self).__init__(arr)
	
		self.fitter = 'sifto'
		self.colour = float('nan')
		self.colour_err = float('nan')
		self.stretch = float('nan')
		self.stretch_err = float('nan')

	#will need to add array details for the SifTO input files
	
	__doc__ = "SiFTO subclass of Supernova Class"
	
#class CosmoChi2:
#	
#	def __init__(self, lc_bmag, lc_bmag_err, lc_width, lc_width_err, lc_colour, lc_colour_err, alpha, beta, omega_m, omega_l, sigma_int):
#		
#		self.lc_bmag = lc_bmag
#		self.lc_bmag_err = lc_bmag_err
#		self.lc_width = lc_width
#		self.lc_width_err = lc_width_err
#		self.lc_colour = lc_colour
#		self.lc_colour_err = lc_colour_err
#		self.alpha = alpha
#		self.beta = beta
#		self.omega_m = omega_m
#		self.omega_l = omega_l
#		self.sigma_int = sigma_int
#		
#	def __call__(self, scriptm):
#		
#		data = 
#


	
def cosmochi2(parameters, zed, bmag=None, bmag_err=None, lc_width=None, lc_width_err=None, lc_colour=None, lc_colour_err=None):
	
	alpha = parameters['alpha'].value
	beta = parameters['beta'].value
	omega_m = parameters['omega_m'].value
	parameters['omega_m'].vary = False
	omega_l = parameters['omega_l'].value
	parameters['omega_l'].vary = False
	scriptm = parameters['scriptm'].value
	int_disp = parameters['int_disp'].value
	
	print omega_m, omega_l
	#pdb.set_trace()
	#measured = data['bmag'] - scriptm + alpha*data['lc_width'] - beta*['lc_colour']
	#measured = bmag - scriptm + alpha*lc_width- beta*lc_colour
	
	
	
	#measured_err = np.sqrt(data['bmag_err']**2  + alpha**2*data['lc_width_err']**2  + beta**2 * data['lc_colour_err'] **2)
	measured_err = np.sqrt(bmag_err**2  + alpha**2*lc_width_err**2  + beta**2 * lc_colour_err **2)

	model = 5.0 * np.log10(np.array([integrate.quad(integrand, 0, zz, args = (omega_m, omega_l))[0] for zz in zed])*1e5)	
	
		#chi2 = sum((data-model)**2./(data_err**2 + int_disp**2))
	chi2 = sum((bmag-(model - alpha*lc_width + beta*lc_colour + scriptm)**2.)/(measured_err**2+ int_disp**2))
	
	return chi2/len(zed)
	

def main():

	parser = argparse.ArgumentParser(description='Simply Python Cosmology Fitter')
	parser.add_argument("filename", type=str, help="Input file of SN data. Name in first column")
	parser.add_argument("style", type=str, choices=['salt2','sifto', 'snoopy' ],help="Style of input data")
	parser.add_argument("sigma", type=float, help='Instrinsic dispersion')
	parser.add_argument("-a", "--alpha", type=float, help='Fix coeff for LC width')
	parser.add_argument("-b", "--beta", type=float, help='Fix coeff for LC colour')
	parser.add_argument("-m", "--omega_m", type=float, help='Fix value of Omega_matter')
	parser.add_argument("-l", "--omega_l", type=float, help='Fix value of Omega_lambda')
	parser.add_argument("-f", "--flat", action = "store_true", help='Force a flat universe model')
	parser.add_argument("-z", "--highz", action="store_true", help='Add high-z data')
	parser.add_argument("-e", "--expansion_rate", type=float, help='Value of the Hubble constant')
	
	args = vars(parser.parse_args())
	
	
	
	
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
			
			if idx == 0:  #first column must be name
				setattr(snvals[dat[0]], col, dat[idx])
			else:
				setattr(snvals[dat[0]], col, float(dat[idx]))
	
	
	#put parameters into class for lmfit
				
	params = Parameters()
	
	try: 
		
		args['alpha']
		
		params.add('alpha', value = args['alpha'], vary=False)
		
	except KeyError:
		
		print 'Alpha free'
		
		params.add('alpha', vary=True, value = 1.0, min=-1, max = 10)
		
	try: 
	
		args['beta']
	
		params.add('beta', value = args['beta'], vary=False)
	
	except KeyError:
	
		print 'Beta free'
	
		params.add('beta', vary=True, value = 1.0, min=-1, max = 10)
		
	try: 
	
		args['omega_m']
	
		params.add('omega_m', value = args['omega_m'], vary=False)
		
		try: 
			args['flat']		
			params.add('omega_l', value = 1.0 - args['omega_m'], vary = False, min=0)
	
		except KeyError:
			print 'Not assuming flat universe'
	
	except KeyError:
	
		print 'Omega_m free'
	
		params.add('omega_m', vary=True, value = 0.25, min=0, max = 5)
		
	try: 
	
		args['omega_l']
	
		params.add('omega_l', value = args['omega_l'], vary=False)
	
		try: 
			args['flat']		
			params.add('omega_m', value = 1.0 - args['omega_l'], vary = False)
	
		except KeyError:
			print 'Not assuming flat universe'
	
	except KeyError:
	
		print 'Omega_l free'
	
		params.add('omega_l', vary=True, value = 0.7, min=0, max = 5)
	
	
	params.add('int_disp', vary = False, value = args['sigma'])
	params.add('scriptm', vary=True, value = 1, min=-70, max = 70)
	
	#Put the data into a dictionary to make it easier
	
	data = dict.fromkeys(['bmag', 'bmag_err', 'zed', 'zed_err', 'lc_width', 'lc_width_err', 'lc_colour', 'lc_colour_err'])
	
	data['bmag'] = np.array([snvals[kk].bmax for kk in snvals.keys()])
	data['bmag_err'] = np.array([snvals[kk].bmax_err for kk in snvals.keys()])
	
	if args['style'] == 'salt2':
		data['lc_width'] = np.array([snvals[kk].x1 for kk in snvals.keys()])
		data['lc_width_err'] = np.array([snvals[kk].x1_err for kk in snvals.keys()])
		
	elif args['style'] =='sifto':
		data['lc_width'] =np.array([snvals[kk].stretch - 1 for kk in snvals.keys()])
		data['lc_width_err'] = np.array([snvals[kk].stretch_err for kk in snvals.keys()])
		
	else:
		print 'SNOOPY STYLE - not yet implemented'
		sys.exit()
		
	data['lc_colour'] = np.array([snvals[kk].colour for kk in snvals.keys()])
	data['lc_colour_err'] = np.array([snvals[kk].colour_err for kk in snvals.keys()])
	
	data['zed'] = np.array([snvals[kk].z for kk in snvals.keys()])
	data['zed_err'] = np.array([snvals[kk].z_err for kk in snvals.keys()])
	
	#pdb.set_trace()
	
	
	
	#result = minimize(cosmochi2, params, args = (data['bmag'], data['bmag_err'], data['lc_width'], data['lc_width_err'], data['lc_colour'], data['lc_colour_err'], data['zed']))
	
	
	bmag = data['bmag']
	bmag_err = data['bmag_err']
	lc_width = data['lc_width']
	lc_width_err = data['lc_width_err']
	lc_colour = data['lc_colour']
	lc_colour_err = data['lc_colour_err']
	zed = data['zed']
	
	
	#result = minimize(cosmochi2, params, args = (bmag, bmag_err, lc_width, lc_width_err, lc_colour, lc_colour_err, zed))
	
	result = minimize(cosmochi2, params, args = (zed,), kws = {'bmag':bmag, 'bmag_err':bmag_err, 'lc_width':lc_width, 'lc_width_err':lc_width_err, 'lc_colour':lc_colour, 'lc_colour_err':lc_colour_err}, method='nelder')
	
	
	report_fit(params)
	
		
	
		
"""
def chi2(some arguments):
	chi2 = sum((data-model)**2./(sigma**2 + int_disp**2))
	
def import_highz():
	work out some way to import the high-z data to add to the minimisation	
"""

if __name__ == "__main__":
	main()