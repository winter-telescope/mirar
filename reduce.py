from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from glob import glob
import os
import subprocess
import warnings
import numpy as np


def reduce(scilist,masterFlatname='masterFlat.fits',masterBiasname='masterBias.fits',filt='r',procDir='proc',calDir='cals',baseDir='data',flatnanThreshold = 0.1):

	if not os.path.exists('%s/%s/%s'%(baseDir,calDir,masterBiasname)):
		print('Bias image not found in ','%s/%s/%s'%(baseDir,calDir,masterBiasname))
		return -1

	img = fits.open('%s/%s/%s'%(baseDir,calDir,masterBiasname))
	masterBias = img[0].data
	img.close()

	proclist = []
	for ind in range(len(scilist)):
		img = fits.open(scilist[ind])
		data = img[0].data
		primaryHeader = img[0].header
		img.close()


		procname = scilist[ind].replace('.fits','.proc.fits')
		procname = '%s/%s/%s'%(baseDir,procDir,procname.split('/')[-1])

		filt = primaryHeader['FILTER']
		if not os.path.exists('%s/%s/%s_%s.fits'%(baseDir,calDir,masterFlatname.replace('.fits',''),filt)):
			print('Flat image not found in ','%s/%s/%s_%s.fits'%(baseDir,calDir,masterFlatname.replace('.fits',''),filt))
			print('Skipping %s'%(scilist[ind]))
			continue

		if primaryHeader['OBSTYPE'] != 'science':
			print('Obstype is not science, skipping %s'%(scilist[ind]))
			continue
		img = fits.open('%s/%s/%s_%s.fits'%(baseDir,calDir,masterFlatname.replace('.fits',''),filt))
		masterFlat = img[0].data
		img.close()

		#print(masterFlat)
		#print(masterBias)
		masterFlatFixed = np.copy(masterFlat)
		if np.any(masterFlat < flatnanThreshold):
			# Set all flat-field values lower than 0.2 to NaN
			masterFlatFixed[masterFlat < flatnanThreshold] = float('NaN')

		'''
		if os.path.exists(procname):
			print('Image is already processed Skipping')
			continue
		'''
		procData = (data - masterBias)/masterFlatFixed
		#print(procData)
		procHDU = fits.PrimaryHDU(procData)  # Create a new HDU with the processed image data
		procHDU.header = primaryHeader       # Copy over the header from the raw file
		procHDU.header.add_history('Bias corrected and flat-fielded') # Add a note to the header
		primaryHeader['BZERO'] = 0
	    # Write the reduced frame to disk
	    
		print(scilist[ind],'->',procname)
		proclist.append(procname)
		procHDU.writeto(procname, overwrite=True)

	return proclist

if __name__== '__main__':
	import argparse
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--d",type=str,default='data',help='Path to the directory where raw images are stored')
	parser.add_argument("--outdir",type=str,default='proc',help='Name of output subdirectory for processed images')
	parser.add_argument("--no_astrom",action="store_false")
	parser.add_argument("--caldir",type=str,default='cals',help='Name of output subdirectory for processed images')
	parser.add_argument("--biasname",type=str,default='masterBias.fits',help='Name of master bias file')
	parser.add_argument("--flatname",type=str,default='masterFlat.fits',help='Name of master flat file')
	parser.add_argument("--filter",type=str,default='u,r',help='Filters, separated by ,')


	args = parser.parse_args()

	imglist = glob('%s/*.fits'%(args.d))

	scilist = []

	filters = args.filter.split(',')
	for filt in filters:
		for imgfile in imglist:
			img = fits.open(imgfile)
			header = img[0].header
			img.close()
	
			if not 'OBSTYPE' in header.keys():
				continue

			if (header['OBSTYPE'] == 'science') & (header['FILTER']==filt):
				scilist.append(imgfile)
	
		print('Found %s sci frames in filter %s'%(len(scilist),filt))
	
		if not os.path.exists('%s/%s'%(args.d,args.outdir)):
			os.mkdir('%s/%s'%(args.d,args.outdir))
			print('Had to make %s directory.'%(args.outdir))

		if not os.path.exists('%s/%s/%s'%(args.d,args.caldir,args.biasname)):
			print('Bias image not found in ','%s/%s/%s'%(args.d,args.caldir,args.biasname))
			exit(0)

		if not os.path.exists('%s/%s/%s_%s.fits'%(args.d,args.caldir,args.flatname.replace('.fits',''),filt)):
			print('Flat image not found in ','%s/%s/%s_%s.fits'%(args.d,args.caldir,args.flatname.replace('.fits',''),filt))
			continue				
	
		reduce(scilist,args.flatname,args.biasname,filt=filt,procDir=args.outdir,calDir=args.caldir,baseDir=args.d)


