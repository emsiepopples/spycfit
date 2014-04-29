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
from lmfit import minimize, Parameters, report_errors
import pdb
from scipy import integrate
from IPython import embed

class Supernova(object):
	
	def __init__(self, arr,lcfitter):
	
		self.name=arr[0]
		self.z = float(arr[1])
		self.z_err = float(arr[2])
		self.maxdate = float(arr[3])
		self.maxdate_err = float(arr[4])
		self.bmax = float(arr[5])
		self.bmax_err = float(arr[6])
		self.type = lcfitter
		
		if self.type == 'salt2':
		
			self.colour = float(arr[11])
			self.colour_err = float(arr[12])
			self.x0 = float(arr[7])
			self.x0_err = float(arr[8])
			self.x1 = float(arr[9])
			self.x1_err = float(arr[10])
			
		elif self.type == 'sifto':
		
	 		self.colour = float('nan')
	 		self.colour_err = float('nan')
	 		self.stretch = float('nan')
	 		self.stretch_err = float('nan')
			
		else:
			
			raise
		
	
		

	
	__doc__ = "Supernova Class"
	

def cosmochisqu(params, snlist):
    #unpack params here, alpha beta etc

    alpha = params['alpha'].value
    beta = params['beta'].value
    omega_m = params['omega_m'].value
    omega_l = params['omega_l'].value
    int_disp = params['int_disp'].value
    scriptm = params['scriptm'].value
    
    #unpack parameters as we did earlier
    redshift = np.array([s.z for s in snlist])
    redshift_err = np.array([s.z_err for s in snlist])
    bmag = np.array([s.bmax for s in snlist])
    bmag_err = np.array([s.bmax_err for s in snlist])
    
    try:
        snlist[0].type == 'salt2'
        lc_width = np.array([s.x1 for s in snlist])
        lc_width_err = np.array([s.x1_err for s in snlist])
        lc_colour = np.array([s.colour for s in snlist])
        lc_colour_err = np.array([s.colour_err for s in snlist])
        
    except:
        snlist[0].type == 'sifto'
        lc_width = np.array([s.stretch - 1.0  for s in snlist])
        lc_width_err = np.array([s.stretch_err for s in snlist])
        lc_colour = np.array([s.colour for s in snlist])
        lc_colour_err = np.array([s.colour_err for s in snlist])
        
    #write this as a separate function in case I want to play with it later
    correction, correction_err = corr_two_params(alpha, beta, [lc_width, lc_width_err], [lc_colour, lc_colour_err])
    #ok this seems to return no errors to here so I think we're ok
    
    model = 5.0 * np.log10(script_lumdist(omega_m, omega_l, redshift))
    
    data = bmag + correction - scriptm
    
    err = np.sqrt(bmag_err**2 + correction_err**2)
    
    chisqu =  (np.sum((model - data)**2) / np.sum(int_disp**2 +err**2))/len(data)
    print chisqu
	
    return (model - data)/np.sqrt((int_disp**2 + err**2))
	

def corr_two_params(aa, bb, width, col):
    correction = aa*width[0] - bb * col[0]
    correction_err = np.sqrt(aa**2 * width[1]**2 + bb**2 * col[1]**2)
    return [correction, correction_err]
	
def integralbit(zz, wm, wl):
    
    return ((1+zz)**2 * (1+wm*zz) - zz*(2+zz)*wl)**-0.5

def script_lumdist(wm, wl, zed):

	zed = np.array(zed)

	if np.abs(1.0 - (wm +wl)) <=0.001:

		curve = 1.0

		const = 299792458. * (1+zed)/np.sqrt(np.abs(curve))
		cosmobit = [np.sqrt(curve) * integrate.quad(integralbit, 0, zz, args=(wm, wl,))[0] for zz in zed]
		return const * cosmobit

	elif wm+wl > 1.:

		print 'Positive curvature'
		curve = 1.0 - wm - wl
		const = 299792458. * (1+zed)/np.sqrt(np.abs(curve))
		cosmobit = np.sin([np.sqrt(np.abs(curve)) * integrate.quad(integralbit, 0, zz, args=(wm,wl,))[0] for zz in zed])
		return const * cosmobit

	else:

		print 'Negative curvature'
		curve = 1.0 - wm - wl
		const = 299792458. * (1+zed)/np.sqrt(np.abs(curve))
		cosmobit = np.sinh([np.sqrt(np.abs(curve)) * integrate.quad(integralbit, 0, zz, args=(wm,wl,))[0] for zz in zed])
		return const * cosmobit

    
 
def three_sigma_clip(snzed, vals, errs, params):

	dist_mod_zed = 5.0 * np.log10(script_lumdist(params['omega_m'].value, params['omega_l'].value, snzed))

	diff = np.abs(vals - dist_mod_zed)/errs  #difference in sigmas

	return np.where(diff >= 3.0)[0]





def main():

	parser = argparse.ArgumentParser(description='Simply Python Cosmology Fitter')
	parser.add_argument("filename", type=str, help="Input file of SN data. Name in first column")
	parser.add_argument("style", type=str, choices=['salt2','sifto', 'snoopy' ],help="Style of input data. NB Snoopy not implemented yet!")
	parser.add_argument("sigma", type=float, help='Instrinsic dispersion')
	parser.add_argument("-a", "--alpha", type=float, help='Fix coeff for LC width')
	parser.add_argument("-b", "--beta", type=float, help='Fix coeff for LC colour')
	parser.add_argument("-m", "--omega_m", type=float, help='Fix value of Omega_matter')
	parser.add_argument("-l", "--omega_l", type=float, help='Fix value of Omega_lambda')
	parser.add_argument("-f", "--flat", action = "store_true", help='Force a flat universe model')
	parser.add_argument("-z", "--highz", action="store_true", help='Add high-z data')
	parser.add_argument("-e", "--expansion_rate", type=float, help='Value of the Hubble constant')
	
	args = vars(parser.parse_args())
	
	print args
	
	
	f = open(args['filename'], 'r')
	lines = f.readlines()
	lines = [x for x in lines if not x.startswith('#')]  #removes comments
	f.close()

	linebyline = [ll.split() for ll in lines] 
	
	sne = [Supernova(lll, args['style']) for lll in linebyline]
	
	#put parameters into class for lmfit
	#this whole section just organises everything from argparse
				
	params = Parameters()
	
	if (args['alpha'] != None):
		
		params.add('alpha', value = args['alpha'], vary=False)
		
	else:
		
		print 'Alpha free'
		
		params.add('alpha', vary=True, value = 1.0, min=-1., max = 10.)
		
	if (args['beta']!= None):
	
		params.add('beta', value = args['beta'], vary=False)
	
	else:
	
		print 'Beta free'
	
		params.add('beta', vary=True, value = 1.0, min=-1, max = 10)
		
	if (args['omega_m']!= None):
	
		params.add('omega_m', value = args['omega_m'], vary=False)
		
		if (args['flat'] == True):

			params.add('omega_l', value = 1.0 - args['omega_m'], vary = False, min=0)
	
		else:
			print 'Not assuming flat universe'
	
	else:
	
		print 'Omega_m free'
	
		params.add('omega_m', vary=True, value = 0.25, min=0, max = 5)
		
	if (args['omega_l']!= None):
	
		params.add('omega_l', value = args['omega_l'], vary=False)
	
		if (args['flat']	!= None	):

			params.add('omega_m', value = 1.0 - args['omega_l'], vary = False)
	
		else:

			print 'Not assuming flat universe'
	
	else:
	
		print 'Omega_l free'
	
		params.add('omega_l', vary=True, value = 0.7, min=0, max = 5)
	
	
	#add these two parameters which are easier

	params.add('int_disp', vary = False, value = args['sigma'])
	params.add('scriptm', vary = True, value = 10.0)
	 

	result = minimize (cosmochisqu, params, args=(sne,))
	cosmochisqu(params, sne)
	report_errors(params)


	test_z = np.arange(0.001, 0.15, 0.001)
	model = 5.0 * np.log10(script_lumdist(params['omega_m'].value, params['omega_l'].value, test_z))
	data = np.array([s.bmax for s in sne]) + params['alpha'].value * np.array([s.x1 for s in sne]) - params['beta'].value * np.array([s.colour for s in sne]) - params['scriptm'].value

	errors = np.sqrt(np.array([s.bmax_err for s in sne])**2 +params['alpha'].value**2 * np.array([s.x1_err for s in sne])**2 \
		+ params['alpha'].stderr**2 * np.array([s.x1 for s in sne]) + params['beta'].value**2 * np.array([s.colour_err for s in sne])**2\
		+ params['beta'].stderr**2 * np.array([s.colour for s in sne])**2 + params['scriptm'].stderr**2)

	

	# units of H0 is wrong. i have no idea what they should be
	#january = np.array([s.bmax for s in sne]) + 0.1 * np.array([s.x1 for s in sne]) - 2.7 * np.array([s.colour for s in sne]) + (-19.08  -5.0 * np.log10(70.*1e3 / 3.0857e22) + 25)

	fig = plt.figure()
	plt.plot(test_z, model, 'k-', label='Best fit')
	#plt.plot(np.array([s.z for s in sne]), data, 'bo', label = 'Data')
	plt.errorbar(np.array([s.z for s in sne]), data, yerr = errors, ecolor='b', fmt='o', mfc='b', label = 'Data')
	#plt.plot(np.array([s.z for s in sne]), january, 'ro', label = 'Jan')
	plt.legend()
	plt.show()
		
	#need to do a 3sigma clipping comparison here

	print 'HERE ARE THE 3SIGMA OUTLIERS'

	clip = three_sigma_clip([s.z for s in sne], data, errors, params)

	for cc in clip: print sne[cc].name



if __name__ == "__main__":
	main()