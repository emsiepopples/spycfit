#!/usr/bin/env python
# encoding: utf-8
"""
weirdness.py

Created by Emma Walker on 2013-10-10.

"""

import sys
import os
import argparse


parser = argparse.ArgumentParser(description='Description of your program')
#parser.add_argument('-f','--foo', help='Description for foo argument', required=True)
#parser.add_argument('-b','--bar', help='Description for bar argument', required=True)
#args = vars(parser.parse_args())
 
#parser = argparse.ArgumentParser()
parser.add_argument("message", type=str, help="Message")
parser.add_argument('-f', '--farg', type=str, help='Test arg1')
parser.add_argument('-d', '--default', type=str, default='DEFAULT_MESSAGE', help='Default test')
args = vars(parser.parse_args())

print args

print args['message']

if args['farg']:
	print 'F arg found'

if args['default']:
	print args['default']
 