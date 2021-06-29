from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.stats import sigma_clip, sigma_clipped_stats
from glob import glob
from astropy.io import fits
import numpy as np
import os
import sys
sys.path.append('.')
from run_astrometry import run_sextractor
import ldactools as aw
import subprocess


astrom_scamp = 'config/scamp.conf'
astrom_sex = 'config/astrom.sex'
astrom_param = 'config/astrom.param'
astrom_filter = 'config/default.conv'
astrom_swarp = 'config/config.swarp'
astrom_nnw = 'config/default.nnw'


def gen_weight_img(imgname):
	img = fits.open(imgname)
	data = img[0].data
	img.close()
	weight_img = np.ones(data.shape)
	weight_img[np.isnan(weight_img)] = 0
	weight_hdu = fits.PrimaryHDU(weight_img)
	weight_hdu.writeto(imgname+'.weight',overwrite=True)


def run_swarp_stack(sciimglist):

	sciimstr = ''
	weightimstr = ''
	for sciimage in sciimglist:
		if not os.path.exists(sciimage+'.weight'):
			gen_weight_img(sciimage)
		sciimstr = sciimstr + sciimage + ','
		weightimstr = weightimstr + sciimage + '.weight' + ','

	sciimstr = sciimstr[:-1]
	weightimstr = weightimstr[:-1]

	img = fits.open(sciimglist[0].replace('.scampastrom',''))
	header = img[0].header
	img.close()

	objname = header['OBJECT']
	filt = header['FILTER']

	w = WCS(header)
	cd11 = header['CD1_1']
	cd21 = header['CD2_1']
	cd12 = header['CD1_2']
	cd22 = header['CD2_2']
	nxpix = header['NAXIS1']
	nypix = header['NAXIS2']
	image_x_cen = nxpix/2
	image_y_cen = nypix/2
	[ra, dec] = w.all_pix2world(image_x_cen,image_y_cen,1)
	xscale = np.sqrt(cd11**2+cd21**2)
	yscale = np.sqrt(cd12**2 + cd22**2)	
	
	pixscale = xscale*3600
	imgpisxize = max(nxpix,nypix)
	stack_outname = '%s_%s_stacked.fits'%(objname,filt)
	#Run swarp on the proc image file
	try:
		#'swarp wasp_ZTF21abcwyvi24.proc.autoastrom.scampastrom.fits,wasp_ZTF21abcwyvi25.proc.autoastrom.scampastrom.fits,wasp_ZTF21abcwyvi29.proc.autoastrom.scampastrom.fits -c config.swarp -RESAMPLE Y -RESAMPLE_DIR . -SUBTRACT_BACK Y -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE wasp_ZTF21abcwyvi24.proc.autoastrom.scampastrom.fits.weight,wasp_ZTF21abcwyvi25.proc.autoastrom.scampastrom.fits.weight,wasp_ZTF21abcwyvi29.proc.autoastrom.scampastrom.fits.weight -COMBINE Y -COMBINE_TYPE MEDIAN -IMAGEOUT_NAME wasp_ZTF21abcwyvi24-25-29.autoastrom.scampastrom.fits -WEIGHTOUT_NAME wasp_ZTF21abcwyvi24-25-29.autoastrom.scampastrom.fits.weight'
		command = 'swarp ' + sciimstr + ' -c ' + astrom_swarp + '  -RESAMPLE Y -RESAMPLE_DIR ' + '.' + ' -SUBTRACT_BACK Y -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE ' + weightimstr + ' ' + ' -COMBINE Y -COMBINE_TYPE MEDIAN -PIXEL_SCALE %.3f'%(pixscale) + ' -IMAGE_SIZE %i'%(imgpisxize) + ' -CENTER %.8f,%.8f'%(ra,dec) + ' -COPY_KEYWORDS FILTER,OBJECT ' + ' -IMAGEOUT_NAME ' + stack_outname +' -WEIGHTOUT_NAME ' + stack_outname.replace('.fits','') + '.resamp.weight -SUBTRACT_BACK N -GAIN 1'
		print('Executing command : %s'%(command))
		rval = subprocess.run(command.split(),check=True,capture_output=True)
		print('Process completed')
		print(rval.stdout.decode())

	except subprocess.CalledProcessError as err:
		print('Could not run swarp with error %s.'%(err))


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--d",type=str,default='data',help='Path to the directory where raw images are stored')
	parser.add_argument("object",type=str,help='Name of object to stack')
	parser.add_argument("--outdir",type=str,default='proc',help='Name of output subdirectory for processed images')
	parser.add_argument("--caldir",type=str,default='cals',help='Name of output subdirectory for processed images')
	parser.add_argument("--biasname",type=str,default='masterBias.fits',help='Name of master bias file')
	parser.add_argument("--flatname",type=str,default='masterFlat.fits',help='Name of master flat file')
	parser.add_argument("filter",type=str,default='r',help='Filter to stack')
	args = parser.parse_args()

	sciimlist = glob('%s/%s/*_%s*.scampastrom.fits'%(args.d,args.outdir,args.object))
	print('Found %s images'%(len(sciimlist)))
	run_swarp_stack(sciimlist)
