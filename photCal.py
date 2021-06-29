from astropy.wcs import WCS
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.stats import sigma_clip, sigma_clipped_stats
from glob import glob
from astropy.io import fits
import numpy as np
import os
import sys
sys.path.append('.')
import subprocess
#from run_astrometry import run_sextractor
import ldactools as aw
import matplotlib.pyplot as plt


photom_sex = 'config/photomCat.sex'
astrom_scamp = 'config/scamp.conf'
astrom_sex = 'config/astrom.sex'
astrom_param = 'config/astrom.param'
photom_param = 'config/photom.param'
astrom_filter = 'config/default.conv'
astrom_swarp = 'config/config.swarp'
astrom_nnw = 'config/default.nnw'

def run_sextractor_phot(imgname,weightimage='weight.fits',gain=1):
	#Run sextractor on the proc image file
	try:
		command = 'sex -c ' + photom_sex + ' ' + imgname + ' ' + '-CATALOG_NAME ' + imgname + '.photom.cat' + ' ' + '-PARAMETERS_NAME ' + photom_param + ' ' + '-FILTER_NAME ' + astrom_filter + ' ' + '-STARNNW_NAME ' + astrom_nnw + ' ' + '-WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE ' + weightimage + ' -CHECKIMAGE_NAME '+imgname+'.bkg,'+imgname+'.bkg.rms -GAIN ' + str(gain)
		print('Executing command : %s'%(command))
		rval = subprocess.run(command.split(),check=True,capture_output=True)
		print('Process completed')
		print(rval.stdout.decode())

	except subprocess.CalledProcessError as err:
		print('Could not run sextractor with error %s.'%(err))



def calc_zeropoint(imgname,boxsize=20,maxmag=18,update_headers=True,weightimage='weight.fits',write_reg=True,ngoodthreshold=5,plot=True):
	img = fits.open(imgname)
	header = img[0].header
	data = img[0].data
	img.close()
	
	w = WCS(header)
	
	[raImage, decImage] = w.all_pix2world(data.shape[0]/2, data.shape[1]/2,1)
	
	catName = 'V/147'
	filt = header['FILTER']
	#filt = 'u'
	print('Querying sdss')
	v = Vizier(columns=['*'], column_filters={"%smag"%(filt):"< %.2f"%(maxmag),"e_%smag"%(filt):"<%.3f"%(1.086/3)}, row_limit=-1)
	Q = v.query_region(SkyCoord(ra=raImage,dec=decImage,unit = (u.deg,u.deg)), radius = str(boxsize)+'m', catalog=catName, cache=False)
	
	if len(Q) == 0:
		print('No matches found in the given radius in SDSS')
		clipped_catalog = []
	else:
		catalog = Q[0]
		print('Retrieved %s sources'%(len(catalog)))
		catalog_crds = [[x['RA_ICRS'],x['DE_ICRS']] for x in catalog]
		with open('sdss.cat.reg','w') as f:
			f.write('wcs\n')
			for c in catalog_crds:
				f.write('CIRCLE(%.4f,%.4f,0.001)\n'%(c[0],c[1]))
		catalog_pixs = w.all_world2pix(catalog_crds,1)
		catalog_pixs = np.array(catalog_pixs)

		catalog_x = catalog_pixs[:,0]
		catalog_y = catalog_pixs[:,1]

		xpix_lolim = 100
		xpix_uplim = 2000
		ypix_lolim = 100
		ypix_uplim = 2000
	
		print(catalog_x,catalog_y)

		edge_mask = (catalog_x<xpix_uplim) & (xpix_lolim<catalog_x) & (ypix_lolim<catalog_y) & (catalog_y<ypix_uplim)
		clipped_catalog = catalog[edge_mask]
		clipped_catalog_crds = SkyCoord(ra=clipped_catalog['RA_ICRS'],dec=clipped_catalog['DE_ICRS'],unit=(u.deg,u.deg))

	print('Found %s good sources in SDSS'%(len(clipped_catalog)))

	if len(clipped_catalog)<ngoodthreshold:
		print('Very few good sources found in SDSS, trying PanStarrs catalog')
		if filt =='u':
			print('Cannot query PanStarrs for u filter :/.')
			return -1
		catName = 'II/349'
		filt = header['FILTER']
		#filt = 'u'
		print('Querying PS1')
		v = Vizier(columns=['*'], column_filters={"%smag"%(filt):"< %.2f"%(maxmag),"e_%smag"%(filt):"<%.3f"%(1.086/3)}, row_limit=-1)
		Q = v.query_region(SkyCoord(ra=raImage,dec=decImage,unit = (u.deg,u.deg)), radius = str(boxsize)+'m', catalog=catName, cache=False)
		catalog = Q[0]
		print(catalog)
		print('Retrieved %s sources'%(len(catalog)))
		catalog_crds = [[x['RAJ2000'],x['DEJ2000']] for x in catalog]
		with open('ps1.cat.reg','w') as f:
			f.write('wcs\n')
			for c in catalog_crds:
				f.write('CIRCLE(%.4f,%.4f,0.001)\n'%(c[0],c[1]))
		catalog_pixs = w.all_world2pix(catalog_crds,1)
		catalog_pixs = np.array(catalog_pixs)
	
		catalog_x = catalog_pixs[:,0]
		catalog_y = catalog_pixs[:,1]
	
		xpix_lolim = 100
		xpix_uplim = 2800
		ypix_lolim = 100
		ypix_uplim = 2800
	
		print(catalog_x,catalog_y)

		edge_mask = (catalog_x<xpix_uplim) & (xpix_lolim<catalog_x) & (ypix_lolim<catalog_y) & (catalog_y<ypix_uplim)
		clipped_catalog = catalog[edge_mask]
		clipped_catalog_crds = SkyCoord(ra=clipped_catalog['RAJ2000'],dec=clipped_catalog['DEJ2000'],unit=(u.deg,u.deg))

		print('Found %s good sources in PanStarrs'%(len(clipped_catalog)))
		if len(clipped_catalog)<ngoodthreshold:
			print('Cannot find enough good sources. Sed.')
			return -1

	if not os.path.exists(imgname+'.photom.cat'):
		run_sextractor_phot(imgname,weightimage=weightimage)

	img_cat = aw.get_table_from_ldac(imgname+'.photom.cat')
	print(img_cat['MAG_APER'][0].shape)
	img_crds = SkyCoord(ra=img_cat['ALPHAWIN_J2000'],dec=img_cat['DELTAWIN_J2000'],unit=(u.deg,u.deg))
	clean_mask = (img_cat['FLAGS']==0) & (img_cat['FWHM_WORLD']<4./3600) & (img_cat['X_IMAGE']>100) & (img_cat['X_IMAGE']<2000) & (img_cat['Y_IMAGE']>100) & (img_cat['Y_IMAGE']<2000)
	clean_img_cat = img_cat[clean_mask]
	clean_img_crds = SkyCoord(ra=clean_img_cat['ALPHAWIN_J2000'],dec=clean_img_cat['DELTAWIN_J2000'],unit=(u.deg,u.deg))

	idx, d2d, d3d = clipped_catalog_crds.match_to_catalog_sky(clean_img_crds)

	match_mask = d2d < 1.0*u.arcsec

	matched_catalog = clipped_catalog[match_mask]
	matched_img = clean_img_cat[idx[match_mask]]

	#apertures = np.array([6.0,10.0,14.0,18.0,22.0 ]) #NEED to change based on requirements
	apertures = np.array([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]) #aperture diameters
	zeropoints = []
	for i in range(len(apertures)):
		offsets = np.ma.array(matched_catalog['%smag'%(filt)] - matched_img['MAG_APER'][:,i])
		cl_offset = sigma_clip(offsets)
		nstars = np.sum(np.invert(cl_offset.mask))
		#print(np.median(cl_offset))
		zp_mean, zp_med, zp_std = sigma_clipped_stats(offsets)
		zeroDict = {'diameter':apertures[i],'zp_mean':zp_mean,'zp_median':zp_med,'zp_std':zp_std,'nstars':nstars,'mag_cat':matched_catalog['%smag'%(filt)][np.invert(cl_offset.mask)],'mag_apers': matched_img['MAG_APER'][:,i][np.invert(cl_offset.mask)]}
		zeropoints.append(zeroDict)


	#mag_auto
	offsets = np.ma.array(matched_catalog['%smag'%(filt)] - matched_img['MAG_AUTO'])
	cl_offset = sigma_clip(offsets,sigma=3)
	nstars = np.sum(np.invert(cl_offset.mask))
	#print(np.median(cl_offset))
	zp_mean, zp_med, zp_std = sigma_clipped_stats(offsets,sigma=1)
	zero_auto = zp_med
	zero_auto_std = zp_std
	zero_auto_nstars = nstars
	zero_auto_mag_cat = matched_catalog['%smag'%(filt)][np.invert(cl_offset.mask)]
	zero_auto_mag_auto = matched_img['MAG_AUTO'][np.invert(cl_offset.mask)]

	if update_headers:
		img = fits.open(imgname,'update')
		for zpvals in zeropoints:
			img[0].header['ZP_%s'%(zpvals['diameter'])] = zpvals['zp_mean']
			img[0].header['ZP_%s_std'%(zpvals['diameter'])] = zpvals['zp_std']
			img[0].header['ZP_%s_nstars'%(zpvals['diameter'])] = zpvals['nstars']

		img[0].header.add_history('Calibrated to SDSS')
		img[0].header['ZP_AUTO'] = zero_auto
		img[0].header['ZP_AUTO_STD'] = zero_auto_std
		img[0].header['ZP_AUTO_nstars'] = zero_auto_nstars

		img[0].header['PHOTCAL'] = 'DONE'
		img.close()

	if plot:
		plt.figure(figsize=(8,8))
		for zpval in zeropoints:
			plt.plot(zpval['mag_cat'],zpval['mag_apers']-zpval['mag_cat'],'.',label=zpval['diameter'])

		plt.plot(zero_auto_mag_cat,zero_auto_mag_auto-zero_auto_mag_cat,'s',label='auto')
		plt.legend()
		plt.savefig(imgname+'_zp.pdf',bbox_inches='tight')

	return 0

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--d",type=str,default='data',help='Path to the directory where raw images are stored')
	parser.add_argument("--procdir",type=str,default='proc',help='Path to the directory where raw images are stored')
	parser.add_argument("--caldir",type=str,default='cals',help='Name of output subdirectory for processed images')
	parser.add_argument("--biasname",type=str,default='masterBias.fits',help='Name of master bias file')
	parser.add_argument("--flatname",type=str,default='masterFlat.fits',help='Name of master flat file')
	parser.add_argument("--single",action="store_true",help='Single image')
	parser.add_argument("--imagename",type=str,help='Single image')
	parser.add_argument("--weight",type=str,default='weight.fits',help='Weight image')

	args = parser.parse_args()

	if args.single:
		calc_zeropoint(args.imagename,weightimage=args.weight)

	imglist = glob('%s/%s/*.autoastrom.scampastrom.resamp.fits'%(args.d,args.procdir))
	for imgname in imglist:
		img = fits.open(imgname)
		header = img[0].header
		img.close()
		if not 'PHOTCAL' in header.keys():
			calc_zeropoint(imgname)