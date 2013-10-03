#/usr/bin/python
"""
Code to amalgamate the standard output of SALT2 into 1 file for analysis in the SPyCoFit tool
Will assume that the results are stored in a data tree by object and subtraction version and that the results files have the default name
'result_salt2.dat'
Output will be 2 files (one for each subtraction method) which can be used to test SPyCoFit.py
"""

import numpy as np
import os
import sys
import argparse
from datetime import datetime as dt

import pdb

class Salt2SN(object):
	
	#initialise object as empty
	
	name=''
	z=float('nan')
	z_err =float('nan')
	maxdate = float('nan')
	maxdate_err = float('nan')
	bmax = float('nan')
	bmax_err = float('nan')
	colour = float('nan')
	colour_err = float('nan')
	x0 = float('nan')
	x0_err = float('nan')
	x1 = float('nan')
	x1_err = float('nan') 
	
	filepath = ''
	
	def getdata(self):
		
		fp = open(self.filepath, 'r')
		
		for i,line in enumerate(fp):
			
			if i==29:
				self.maxdate = line.split()[1]
				self.maxdate_err = line.split()[2]
			if i==30:
				self.z = line.split()[1]
				tmperr = line.split()[2]
				if tmperr == '0':
					self.z_err = '0.001'
				else:
					self.z_err = tmperr
			if i==31:
				self.colour = line.split()[1]
				self.colour_err = line.split()[2]
			if i==32:
				self.x0 = line.split()[1]
				self.x0_err = line.split()[2]
			if i==33:
				self.x1 = line.split()[1]
				self.x1_err = line.split()[2]
			if i==34:
				self.bmax = line.split()[1]	
				self.bmax_err = line.split()[2]
		return self
		
	def writeout(self, outfile):
		
		outfile.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(self.name, self.z, self.z_err, self.maxdate, self.maxdate_err, self.bmax, self.bmax_err, self.x0, self.x0_err, self.x1, self.x1_err, self.colour, self.colour_err))
		return
		
		
		
def main(argv=None):
	
	parser = argparse.ArgumentParser()
	parser.add_argument("startdir", type=str, help="Starting directory")
	parser.add_argument("version", type=str, help="Subtraction Version")
	parser.add_argument("-o", "--outfile", type=str, help="Output file", default = 'salt2.dat')
	args = parser.parse_args()
	
	saltname = 'result_salt2.dat'
	
	#As this is for SALT2 we can define the names of all the column headers
	
	col_names = ['name', 'z', 'z_err', 'maxdate', 'maxdate_err', 'bmax', 'bmax_err', 'colour', 'colour_err', 'x0', 'x0_err', 'x1', 'x1_err']
	
	ncol = len(col_names)
	
	fmt = ['S12', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8']
	
	datarr = np.zeros(ncol, dtype={'names':col_names, 'formats':fmt})
	
	#open output file
	g = open(args.outfile, 'w')
	
	
	g.write("\t".join(col_names))
	g.write("\n")
	
	#list objects in the starting directory
	
	sne=os.listdir(args.startdir)
	
	for s in sne:
		
		locate = os.path.join(args.startdir, s, args.version, saltname)
		
		if os.path.isfile(locate):
			sn = Salt2SN()
			sn.filepath = locate
			sn.name = s
			sn.getdata() #use the class method to read date in from file
			sn.writeout(g) #have to give it the file object
			
			
	now = dt.now()
	timestamp = now.strftime("%Y-%m-%d %H:%M")
	
	g.write("# Created on {0}".format(timestamp))
	
	g.close()
	
	
	
	
if __name__=="__main__":
	sys.exit(main())
