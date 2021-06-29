from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from glob import glob
import os
import subprocess
import warnings
import numpy as np


def makeMasterBias(biaslist, xlolim=200,xuplim=1300,ylolim=200,yuplim=1300,calDir='cals',basedir='data',outname='masterBias.fits'):

	img = fits.open(biaslist[0])
	header = img[0].header
	img.close()
	nx = header['NAXIS1']
	ny = header['NAXIS2']

	biases = np.zeros((ny,nx,len(biaslist)))
	for ind in range(len(biaslist)):
		print('Reading bias %i/%i'%(ind+1,len(biaslist)))
		img = fits.open(biaslist[ind])
		biases[:,:,ind] = img[0].data
		img.close()

	print('Median combining %s biases'%(len(biaslist)))

	masterBias = np.nanmedian(biases,axis=2)

	img = fits.open(biaslist[0])
	primaryHeader = img[0].header
	img.close()
	procHDU = fits.PrimaryHDU(masterBias)  # Create a new HDU with the processed image data
	procHDU.header = primaryHeader       # Copy over the header from the raw file
	procHDU.header.add_history('median stacked bias')
	#procHDU.writeto('/media/viraj/New Volume/winterdrp/wasp_data/masterBias.fits',overwrite=True)
	procHDU.writeto('%s/%s/%s'%(basedir,calDir,outname),overwrite=True)
	return 0


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--d",type=str,default='data',help='Path to the directory where flats are stored')
	parser.add_argument("--outdir",type=str,default='cals',help='Name of output subdirectory')
	parser.add_argument("--redo",action="store_true",help='Redo bias')
	
	args = parser.parse_args()

	imglist = glob('%s/*.fits'%(args.d))

	biaslist = []
	for imgfile in imglist:
		img = fits.open(imgfile)
	
		if img[0].header['OBJECT'] == 'bias' :
			biaslist.append(imgfile)

	print('Found %s bias frames'%(len(biaslist)))

	if not os.path.exists('%s/%s'%(args.d,args.outdir)):
		os.mkdir('%s/%s'%(args.d,args.outdir))
		print('Had to make cals directory.')


	if os.path.exists('%s/%s/masterBias.fits'%(args.d,args.outdir)) and not args.redo:
		print('Bias already exists. Rerun with --redo option')
		exit(0)

	makeMasterBias(biaslist,calDir=args.outdir,basedir=args.d)
	print('Wrote master bias.')