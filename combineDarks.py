from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from glob import glob
import os
import subprocess
import warnings
import numpy as np


def makeMasterDarks(darklist, masterbiasname = 'masterBias.fits', xlolim=500,xuplim=3500,ylolim=500,yuplim=3500,calDir='cals',basedir='data',outbasename='masterDark',normalise=False):

	img = fits.open(darklist[0])
	header = img[0].header
	img.close()
	nx = header['NAXIS1']
	ny = header['NAXIS2']
	exptime = header['EXPTIME']

	darks = np.zeros((ny,nx,len(darklist)))

	img = fits.open('%s/%s/%s'%(basedir,calDir,masterbiasname))
	masterBias = img[0].data
	header = img[0].header
	img.close() 

	if not normalise:
		for ind in range(len(darklist)):
			print('Reading dark %i/%i'%(ind+1,len(darklist)))
			img = fits.open(darklist[ind])
			if img['EXPTIME'] != exptime:
				print('Error, all dark files do not have the same exposure time. Please run with normalise option to create a normalised master dark.')
				return -1

			darks[:,:,ind] = img[0].data - masterBias
			img.close()

		print('Median combining %s darks'%(len(darks)))
		masterDarkname = '%s_%s.fits'%(outbasename,exptime)

	if normalise:
		for ind in range(len(darklist)):
			print('Reading dark %i/%i'%(ind+1,len(darklist)))
			img = fits.open(darklist[ind])
			dark_exptime = img[0].header['EXPTIME']
			
			darks[:,:,ind] = (img[0].data - masterBias)/dark_exptime
			img.close()

		print('Median combining %s darks'%(len(darks)))
		masterDarkname = '%s_normalised.fits'%(outbasename)

	masterDark = np.nanmedian(darks,axis=2)


	img = fits.open(darklist[0])
	primaryHeader = img[0].header
	img.close()
	procHDU = fits.PrimaryHDU(masterBias)  # Create a new HDU with the processed image data
	procHDU.header = primaryHeader       # Copy over the header from the raw file
	if normalise:
		procHDU.header.add_history('median stacked normalised dark')
		procHDU.header['EXPTIME'] = 1.0

	else:
		procHDU.header.add_history('median stacked dark')
		procHDU.header['EXPTIME'] = exptime
	#procHDU.writeto('/media/viraj/New Volume/winterdrp/wasp_data/masterBias.fits',overwrite=True)
	procHDU.writeto('%s/%s/%s'%(basedir,calDir,masterDarkname),overwrite=True)
	return 0