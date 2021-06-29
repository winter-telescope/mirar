from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from glob import glob
import os
import subprocess
import warnings
import numpy as np


def makeMasterFlat(flatlist, masterbiasname,xlolim=200,xuplim=1300,ylolim=200,yuplim=1300,calDir='cals',basedir='data'):
	#print('Flatlist',flatlist)
	img = fits.open('%s/%s/%s'%(basedir,calDir,masterbiasname))
	masterBias = img[0].data
	header = img[0].header
	img.close()


	nx = header['NAXIS1']
	ny = header['NAXIS2']

	flats = np.zeros((ny,nx,len(flatlist)))
	print('Combining flats')
	for ind in range(len(flatlist)):
		print('Reading flat %i/%i'%(ind+1,len(flatlist)))
		img = fits.open(flatlist[ind])
		data = img[0].data
		header = img[0].header
		img.close()
		data = data - masterBias
		median = np.nanmedian(data[xlolim:xuplim,ylolim:yuplim])
		flats[:,:,ind] = data/median
		print(flatlist[ind],median,header['EXPTIME'],header['FILTER'])

	print('Median combining %s flats'%(len(flatlist)))

	masterFlat = np.nanmedian(flats,axis=2)
	
	img = fits.open(flatlist[0])
	primaryHeader = img[0].header
	img.close()
	procHDU = fits.PrimaryHDU(masterFlat)  # Create a new HDU with the processed image data
	procHDU.header = primaryHeader       # Copy over the header from the raw file
	procHDU.header.add_history('Stacked flat-fielded')

	filt = primaryHeader['FILTER']
	primaryHeader['BZERO'] = 0
	print('Writing to','%s/%s/masterFlat_%s.fits'%(basedir,calDir,filt))
	procHDU.writeto('%s/%s/masterFlat_%s.fits'%(basedir,calDir,filt),overwrite=True)
	return 0


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--filter",type=str,default='u,r',help='Filters, separated by ,')
	parser.add_argument("--d",type=str,default='data',help='Path to the directory where flats are stored')
	parser.add_argument("--outdir",type=str,default='cals',help='Name of output subdirectory')
	parser.add_argument("--biasname",type=str,default='masterBias.fits',help='Name of master bias file')
	parser.add_argument("--redo",action="store_true",help='Redo flats')
	
	args = parser.parse_args()
	filterlist = args.filter.split(',')

	imglist = glob('%s/*.fits'%(args.d))
	
	for filt in filterlist:
		print('Doing filter %s'%(filt))
		if os.path.exists('%s/%s/masterFlat_%s.fits'%(args.d,args.outdir,filt)) and not args.redo:
			print('Flat already exists, run with --redo option')
			continue

		flatlist = []
		for imgfile in imglist:
			img = fits.open(imgfile)
		
			if (img[0].header['OBJECT'] == 'flat') & (img[0].header['FILTER']==filt):
				flatlist.append(imgfile)

		print('Found %s flats for filter %s'%(len(flatlist),filt))

		if not os.path.exists('%s/%s'%(args.d,args.outdir)):
			os.mkdir('%s/%s'%(args.d,args.outdir))
			print('Had to make cals directory. Please make a masterBias and put it in this directory.')
			exit(0)

		if not os.path.exists('%s/%s/%s'%(args.d,args.outdir,args.biasname)):
			print('Bias file not found, please generate it first.')
			exit(0)


		if (len(flatlist)>0):
			makeMasterFlat(flatlist,args.biasname,calDir=args.outdir,basedir=args.d)
			print('Wrote master flat for filter %s'%(filt))

		else:
			print('No flats found for filter %s'%(filt))