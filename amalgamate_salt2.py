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

def getdata(filepath):
	
	#use this to amalgamate a set of SALT2 files
	
	fp = open(filepath, 'r')
	
	for i,line in enumerate(fp):
		
		if i==29:
			maxdate = line.split()[1]
			maxdate_err = line.split()[2]
		if i==30:
			z = line.split()[1]
			tmperr = line.split()[2]
			if tmperr == '0':
				z_err = '0.001'
			else:
				z_err = tmperr
		if i==31:
			colour = line.split()[1]
			colour_err = line.split()[2]
		if i==32:
			x0 = line.split()[1]
			x0_err = line.split()[2]
		if i==33:
			x1 = line.split()[1]
			x1_err = line.split()[2]
		if i==34:
			bmax = line.split()[1]	
			bmax_err = line.split()[2]
			
	print [[z, z_err], [maxdate, maxdate_err], [bmax, bmax_err],[x0, x0_err], [x1, x1_err], [colour, colour_err]]
	return [[z, z_err], [maxdate, maxdate_err], [bmax, bmax_err],[x0, x0_err], [x1, x1_err], [colour, colour_err]]
	

	
	
		
def main(argv=None):
	
	parser = argparse.ArgumentParser(description='Stand-alone programme to organise the LSQ SALT2 directories and results into 1 input file')
	parser.add_argument("startdir", type=str, help="Starting Directory")
	parser.add_argument('version', type=str, help='Version of subtraction')
	parser.add_argument('-o', '--outfile', type=str, default='salt2_amalgamated.dat', help='Output File')
	args = vars(parser.parse_args())
	
	print args
	
	saltname = 'result_salt2.dat'
	
	#As this is for SALT2 we can define the names of all the column headers
	
	col_names = ['name', 'z', 'z_err', 'maxdate', 'maxdate_err', 'bmax', 'bmax_err', 'x0', 'x0_err', 'x1', 'x1_err', 'colour', 'colour_err']
	
	ncol = len(col_names)
	
	
	#open output file
	g = open(args['outfile'], 'w')
	
	
	g.write("\t".join(col_names))
	g.write("\n")
	
	#list objects in the starting directory
	
	sne=os.listdir(args['startdir'])
	
	
	for s in sne:
		locate = os.path.join(args['startdir'], s, args['version'], saltname)
		
		if os.path.isfile(locate):
	
			name = s
			vals = getdata(locate) 
			g.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(name, vals[0][0], vals[0][1], vals[1][0], vals[1][1], vals[2][0], vals[2][1], vals[3][0], vals[3][1], vals[4][0], vals[4][1], vals[5][0], vals[5][1]))
	
			
	now = dt.now()
	timestamp = now.strftime("%Y-%m-%d %H:%M")
	
	g.write("# Created on {0}".format(timestamp))
	
	g.close()
		
		
		
	
if __name__=="__main__":
	sys.exit(main())
