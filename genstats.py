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


astrom_scamp = 'config/scamp.conf'
astrom_sex = 'config/astrom.sex'
astrom_param = 'config/astrom.param'
astrom_filter = 'config/default.conv'
astrom_swarp = 'config/config.swarp'
astrom_nnw = 'config/default.nnw'
photom_sex = 'config/photomCat.sex'


def run_sextractor(imgname,pixscale=0.47,regions=True,weightimg='weight.fits'):
	#Run sextractor on the proc image file
	try:
		command = 'sex -c ' + astrom_sex + ' ' + imgname + ' ' + '-CATALOG_NAME ' + imgname + '.cat' + ' -CATALOG_TYPE FITS_LDAC ' + '-PARAMETERS_NAME ' + astrom_param + ' ' + '-FILTER_NAME ' + astrom_filter + ' ' + '-STARNNW_NAME ' + astrom_nnw + ' ' + '-WEIGHT_TYPE NONE -CHECKIMAGE_TYPE NONE -PIXEL_SCALE ' + str(pixscale) + ' -DETECT_THRESH 10 -ANALYSIS_THRESH 10 -SATUR_LEVEL 60000 -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE ' + weightimg + ' -CHECKIMAGE_NAME '+imgname+'.seg,'+imgname+'.bkg,'+imgname+'.bkg.rms'
		print('Executing command : %s'%(command))
		rval = subprocess.run(command.split(),check=True,capture_output=True)
		print('Process completed')
		print(rval.stdout.decode())

	except subprocess.CalledProcessError as err:
		print('Could not run sextractor with error %s.'%(err))
		return

	if regions:
		t = aw.get_table_from_ldac(imgname+'.cat')
		with open(imgname+'.cat'+'.stats.reg','w') as f:
			f.write('image\n')
			for row in t:
				f.write('CIRCLE(%s,%s,%s) # text={%.2f,%.2f}\n'%(row['X_IMAGE'],row['Y_IMAGE'],row['FWHM_IMAGE']/2,row['FWHM_IMAGE'],row['SNR_WIN']))




def get_img_fwhm(imgname,pixscale=0.47,weightimg='weight.fits',xlolim=10,xuplim=2000,ylolim=10,yuplim=2000,exclude=False):
	if not os.path.exists(imgname+'.cat'):
		run_sextractor(imgname,pixscale,weightimg=weightimg)

	img_cat = aw.get_table_from_ldac(imgname+'.cat')
	print('Found %s sources'%(len(img_cat)))

	center_mask = (img_cat['X_IMAGE']<xuplim) & (img_cat['X_IMAGE']>xlolim) & (img_cat['Y_IMAGE']<yuplim) & (img_cat['Y_IMAGE']>ylolim)
	if exclude :
		center_mask = np.invert(center_mask)
	print('Using %s sources'%(len(img_cat[center_mask])))
	mean, median, std = sigma_clipped_stats(img_cat[center_mask]['FWHM_IMAGE'])
	return mean,median, std


def gen_map(imgname,pixscale=0.466,weightimg='weight.fits',regions=False):
	npix = 2000
	x_inds = np.linspace(0,npix,5)
	y_inds = np.linspace(0,npix,5)

	if not os.path.exists(imgname+'.cat'):
		run_sextractor(imgname)
	t = aw.get_table_from_ldac(imgname+'.cat')

	fwhms = np.zeros((4,4))
	elongs = np.zeros((4,4))
	for i in range(len(x_inds)-1):
		for j in range(len(y_inds)-1):
			cut_t = (x_inds[i]<t['X_IMAGE'])& (t['X_IMAGE']<x_inds[i+1]) & (y_inds[j]<t['Y_IMAGE'])& (t['Y_IMAGE']<y_inds[j+1])
			mean, median, std = sigma_clipped_stats(t[cut_t]['FWHM_IMAGE'])
			fwhms[j][i] = median
			mean, median, std = sigma_clipped_stats(t[cut_t]['ELONGATION'])
			elongs[j][i] = median


	fig = plt.figure()
	gs = GridSpec(1,2,wspace=0.5)
	ax1 = fig.add_subplot(gs[0])
	fmap = ax1.imshow(fwhms*pixscale,origin='lower',extent=(0,npix,0,npix),vmin=1.8,vmax=3.5)
	ax1.set_xticks(x_inds)
	ax1.set_yticks(y_inds)

	divider = make_axes_locatable(ax1)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	fig.colorbar(fmap, cax=cax)
	ax1.set_title('FWHM')

	ax2 = fig.add_subplot(gs[1])
	emap = ax2.imshow(1-1/elongs,origin='lower',extent=(0,npix,0,npix))
	ax2.set_xticks(x_inds)
	ax2.set_yticks(y_inds)
	divider = make_axes_locatable(ax2)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	ax2.set_title('Ellipticity')
	fig.colorbar(emap, cax=cax)

	plt.savefig('%s_fwhm_ellipticity_maps.pdf'%(imgname.split('/')[-1]),bbox_inches='tight')


def gen_weight_img(imgname):
	img = fits.open(imgname)
	data = img[0].data
	img.close()
	weight_img = np.ones(data.shape)
	weight_img[np.isnan(data)] = 0
	#weight_img[2048,:] = 0	 #Better way to flag this out 
	weight_hdu = fits.PrimaryHDU(weight_img)
	weight_hdu.writeto(imgname+'.weight',overwrite=True)


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--imgname",type=str,help='Imagename')
	parser.add_argument("--pixscale",type=float,default=0.47,help='Pixelscale')
	parser.add_argument("--plot",action='store_true',help='Plot?')
	parser.add_argument("--weight",type=str,default='weight.fits',help='weight image')
	
	args = parser.parse_args()

	mean_fwhm, med_fwhm, std = get_img_fwhm(args.imgname,pixscale=args.pixscale,weightimg=args.weight)
	print('Mean, Median, Std : %.2f,%.2f,%.2f'%(mean_fwhm*args.pixscale,med_fwhm*args.pixscale,std*args.pixscale))

	if args.plot:
		gen_map(args.imgname,pixscale=args.pixscale)

