from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from glob import glob
import os
import subprocess
import warnings
import numpy as np
import ldactools as aw
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import argparse
from reduce import reduce
from combineBias import makeMasterBias
from combineFlats import makeMasterFlat
from run_astrometry import run_astrometry



def prepare(basedir,redo=False,caldir='cals',procdir='proc',otherdir='other',masterBiasname='masterBias.fits',masterFlatname='masterFlat'):
	#Go through all files in a directory. First, check if there is a cals subdirectory with a masterbias and a masterflat for both filters.
	if not os.path.exists(os.path.join(basedir,caldir)):
		print('Cals directory not found, making it now.')
		os.mkdir(os.path.join(basedir,caldir))

	if not os.path.exists(os.path.join(basedir,procdir)):
		print('proc directory not found, making it now.')
		os.mkdir(os.path.join(basedir,procdir))

	if not os.path.exists(os.path.join(basedir,otherdir)):
		print('proc directory not found, making it now.')
		os.mkdir(os.path.join(basedir,otherdir))

	imglist = glob('%s/*.fits'%(basedir))
	print('Found %s images.'%(len(imglist)))
	biaslist = []
	flatlist = {} 
	flatlist['u'] = []
	flatlist['r'] = []
	scilist = []
	sci_filters = []

	filterlist = ['u','r']

	if os.path.exists('%s/%s/%s'%(basedir,caldir,masterBiasname)):
		print('Found masterbias. Will use this.')
		biasexists = 1

	else:
		print('Did not find master bias. Looking for biaslist.txt file')
		if os.path.exists('%s/biaslist.txt'%(basedir)) and not redo:
			biaslist = np.genfromtxt('%s/biaslist.txt'%(basedir))
			if len(biaslist)>0:
				biasexists = 1
			else:
				biasexists = 0
		else:
			print('Did not find biaslist.txt file, or redo option is True. Will search for raw bias files again')
			biasexists = 0


	for filt in filterlist:
		if not os.path.exists('%s/%s/%s_%s.fits'%(basedir,caldir,masterFlatname,filt)):
			print('%s band masterflat not found. Will search for raw flat images in this band'%(filt))
			flatexists = 0

			if os.path.exists('%s/flatlist_%s.txt'%(basedir,filt)) and not redo:
				flatlist[filt] = np.genfromtxt('%s/flatlist_%s.txt'%(basedir,filt))
				if len(flatlist[filt])>0:
					flatexists = 1
				else:
					flatexists = 0

			else:
				print('Did not find flatlist_%s.txt or redo is set on. Will search for flat files again.'%(filt))

		else:
			print('%s band masterflat found! Will use this flat.')
			flatexists = 1

	if biasexists==0 or flatexists==0:
		biaslist = []
		flatlist = {} 
		flatlist['u'] = []
		flatlist['r'] = []
		print('Searching the entire folder for bias and flats. Will also look for science files while at it :)')

		for imgname in imglist:
			img = fits.open(imgname)
			header = img[0].header
			img.close()

			#print(flatlist[header['FILTER']])
			if header['OBSTYPE'] == 'flat':
				flatlist[header['FILTER']].append(imgname)

			elif header['OBSTYPE'] == 'bias':
				biaslist.append(imgname)

			elif header['OBSTYPE'] == 'science':
				scilist.append(imgname)
				sci_filters.append(header['FILTER'])


			else:
				print('Image %s is not a bias, science or a flat. Moving to other directory.'%(imgname))
				os.system('mv %s %s/%s/'%(imgname,basedir,otherdir))

		print('Found %s bias frames'%(len(biaslist)))
		with open('%s/biaslist.txt'%(basedir),'w') as f:
			np.savetxt(f,biaslist,fmt='%s')

		if len(biaslist)==0:
			print('No bias frames found in raw')
			if biasexists==0:
				print('No masterbias found either. Please include bias frames or add a masterBias.fits to the cals folder')	
				return -1
		else:
			print('Combining bias frames')
			makeMasterBias(biaslist,basedir=basedir)


		for filt in filterlist:
			print('Found %s flat frames in %s band'%(len(flatlist[filt]),filt))
			#print(flatlist[filt])
			with open('%s/flatlist_%s.txt'%(basedir,filt),'w') as f:
				np.savetxt(f,flatlist[filt],fmt='%s')
			if len(flatlist[filt])==0:
				print('No flats found for band %s. Will skip reducing observations in this band'%(filt))
			else:
				print('Combining flats for band %s'%(filt))
				makeMasterFlat(flatlist[filt],masterbiasname='masterBias.fits',basedir=basedir)


		print('Found %s science frames'%(len(scilist)))
		print('Science observations are in the following filters',np.unique(sci_filters))

	else:
		allflatlist = []
		for filt in filterlist:
			allflatlist.append(flatlist[filt])
		allflatlist = np.array(allflatlist)
		nonscilist = np.union1d(biaslist,allflatlist.flatten())
		scilist = np.setdiff1d(imglist,nonscilist)

	if not redo:
		proclist = glob('%s/%s/*proc.autoastrom.scampastrom.fits'%(basedir,procdir))
		proclist = [x.replace('.proc.autoastrom.scampastrom','') for x in proclist]
		proclist = [x.replace('/proc','') for x in proclist]
		
		print('Already processed %s files'%(len(proclist)))
		scilist = np.setdiff1d(scilist,proclist)

	print('Will reduce %s science frames'%(len(scilist)))

	proclist = reduce(scilist,baseDir=basedir)

	 
	run_astrometry(proclist)

	
def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--d",type=str,help='Dirname')

	parser.add_argument("--pixscale",type=float,default=0.466,help='Pixelscale')	
	parser.add_argument("--exclude",action='store_true',help='Exclude central values')
	parser.add_argument("--redo",action='store_true')
	args = parser.parse_args()

	imglist = glob('%s/*.fits'%(args.d))
	prepare(args.d,redo=args.redo)

if __name__ == '__main__':
	main()