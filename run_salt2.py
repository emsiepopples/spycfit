#/usr/bin/python
"""Python script to run SALT2 on a set of subdirectories for each object"""

import os
import sys
import glob
import argparse

def main(argv=None):

	parser = argparse.ArgumentParser(description='Runs SALT2 on a set of subdirectories')
	parser.add_argument("startdir", type=str, help="Starting Directory")
	args = vars(parser.parse_args())

	startlocation = os.getcwd()

	topdir =args['startdir']

	for root, dirs, files in os.walk(topdir):

		os.chdir(root)

		try:

			os.rename('result_salt2.dat', 'result_salt22.old')

		except:

			pass

		photfiles = ' '.join(glob.glob("LSQ*.dat"))

		os.system('snfit -w 3000 10000 {0}'.format(photfiles))


	os.chdir(startlocation)


if __name__=="__main__":
	sys.exit(main())