#!/usr/bin/env python
# encoding: utf-8
"""
moveimagefiles.py

this python script moves subset of images from the larger Human_vs_Machine directory of images.
the directory has 8 categories of object, 8 exemplars per category, with 7 variation levels each.

Created by labuser on 2010-08-31.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import glob
import shutil

def main():
	if len(sys.argv) < 2:
		print "moveimagefile.py <dirfrom> <dirto>"
		return
	dirname = sys.argv[1]
	dirto = sys.argv[2]
	objdirs =  glob.glob(dirname+'/*/')
	print objdirs
	flind = []
	#flind has indices the pick 10 our of every 40 images per object
	for ind in range(0,32,4):
		flind = flind+range(ind*10,(ind+1)*10)
	for objdir in objdirs:
		fname = glob.glob(objdir+'*.png')
		catname = objdir.split('/')[-2]
		if 'Variation00' in dirname:
			[shutil.copy(fname[ind],dirto+'/'+catname+'_'+ os.path.basename(fname[ind])) for ind in flind]
		else:
			[shutil.copy(fn,dirto+'/'+catname+'_'+os.path.basename(fn)) for fn in fname]
	print dirname
	print dirto
	pass


if __name__ == '__main__':
	main()

