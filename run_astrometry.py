from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from glob import glob
import os
import subprocess
import warnings
import numpy as np
from astroquery.gaia import Gaia
import sys
sys.path.append('.')
import ldactools as aw
from astropy.io.fits import Header


astrom_scamp = 'config/scamp.conf'
astrom_sex = 'config/astrom.sex'
astrom_param = 'config/astrom.param'
astrom_filter = 'config/default.conv'
astrom_swarp = 'config/config.swarp'
astrom_nnw = 'config/default.nnw'


def secondPass_astrometry(imgname,ra=None,dec=None,catalog_min_mag=10,catalog_max_mag=18,catalog_box_size_arcmin=30,writefiles=True,weightimg='weight.fits'):
	#Run scamp to generate a refined solution
	res = run_sextractor(imgname,weightimg=weightimg)
	if res<0:
		return -1
			
	if ra is None and dec is None:
		img = fits.open(imgname)
		header = img[0].header
		data = img[0].data
		img.close()

		w = WCS(header)
		[ra, dec] = w.all_pix2world(data.shape[0]/2, data.shape[1]/2,1)

	print(ra,dec)
	res = make_gaia_catalog(ra,dec,imgname.replace('.fits','.gaia'),catalog_box_size_arcmin,catalog_min_mag,catalog_max_mag,writeldac=True)
	if res<0:
		return -1

	print('Downloaded gaia sources')
	
	res = run_scamp(imgname)
	if res<0:
		return -1

	#print('Scamp run successfully')
	if writefiles:
		headername = imgname + '.head'
		write_scamped_file(imgname,headername)

	return 0


def write_scamped_file(imgname,headername):
	img = fits.open(imgname)
	header = img[0].header
	data = img[0].data
	img.close()

	with open(headername,'r') as f:
		a = f.read()

	h = Header()
	h = h.fromstring(a,sep='\n')
	for k in h.keys():
		if k=='HISTORY' or k=='COMMENT':
			continue

		header[k] = h[k]
	header.add_history('%s'%(h['HISTORY']))
	procHDU = fits.PrimaryHDU(data)
	procHDU.header = header
	procHDU.writeto(imgname.replace('.fits','')+'.scampastrom.fits',overwrite=True)


def make_gaia_catalog(ra, dec, tmcatname, catalog_box_size_arcmin, catalog_min_mag, catalog_max_mag, writeldac = True):
    print("SELECT * FROM gaiadr2.gaia_source AS g, gaiadr2.tmass_best_neighbour AS tbest, gaiadr1.tmass_original_valid AS tmass WHERE g.source_id = tbest.source_id AND tbest.tmass_oid = tmass.tmass_oid AND CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', %.4f, %.4f, %.4f))=1 AND tmass.j_m > %.2f AND tmass.j_m < %.2f AND tbest.number_of_mates=0 AND tbest.number_of_neighbours=1;"%(ra, dec, catalog_box_size_arcmin/60, catalog_min_mag, catalog_max_mag))
    try:
        job = Gaia.launch_job_async("SELECT * FROM gaiadr2.gaia_source AS g, gaiadr2.tmass_best_neighbour AS tbest, gaiadr1.tmass_original_valid AS tmass WHERE g.source_id = tbest.source_id AND tbest.tmass_oid = tmass.tmass_oid AND CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', %.4f, %.4f, %.4f))=1 AND tmass.j_m > %.2f AND tmass.j_m < %.2f AND tbest.number_of_mates=0 AND tbest.number_of_neighbours=1;"%(ra, dec, catalog_box_size_arcmin/60, catalog_min_mag, catalog_max_mag), dump_to_file = False)
        print('Yay')
        t = job.get_results()
        
    except:
        print('Could not query Gaia .. ')
        return -1
   
    print(t.colnames)
    #t.remove_columns(['designation', 'phot_variable_flag', 'datalink_url', 'epoch_photometry_url', 'original_ext_source_id', 'designation_2'])
    t.remove_columns(['designation', 'phot_variable_flag', 'datalink_url', 'original_ext_source_id', 'designation_2'])
    t['ph_qual'] = t['ph_qual'].astype(str)
    t['ra_errdeg'] = t['ra_error'] / 3.6e6
    t['dec_errdeg'] = t['dec_error'] / 3.6e6
    t['FLAGS'] = 0
       
    if writeldac:
        if os.path.exists(tmcatname + '.ldac'):
            os.remove(tmcatname + '.ldac')
        aw.save_table_as_ldac(t, tmcatname + '.ldac')

    return 0

def run_sextractor(imgname,weightimg='weight.fits'):
	#Run sextractor on the proc image file
	try:
		command = 'sex -c ' + astrom_sex + ' ' + imgname + ' ' + '-CATALOG_NAME ' + imgname + '.cat' + ' ' + '-PARAMETERS_NAME ' + astrom_param + ' ' + '-FILTER_NAME ' + astrom_filter + ' ' + '-STARNNW_NAME ' + astrom_nnw + ' ' + '-WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+ weightimg +' -CHECKIMAGE_TYPE NONE'
		print('Executing command : %s'%(command))
		rval = subprocess.run(command.split(),check=True,capture_output=True)
		print('Process completed')
		print(rval.stdout.decode())
		return 0

	except subprocess.CalledProcessError as err:
		print('Could not run sextractor with error %s.'%(err))
		return -1


def run_scamp(imgname):
	#Run scamp on the proc image file
	try:
		command = 'scamp -c ' + astrom_scamp + ' ' + imgname + '.cat' + ' -ASTREFCAT_NAME '  + imgname.replace('.fits','.gaia') +'.ldac' + ' -ASTREFMAG_KEY phot_rp_mean_mag'
		print('Executing command : %s'%(command))
		rval = subprocess.run(command.split(),check=True,capture_output=True)
		print('Process completed')
		print(rval.stdout.decode())
		return 0

	except subprocess.CalledProcessError as err:
		print('Could not run scamp with error %s.'%(err))
		return -1


def gen_weight_img(imgname):
	img = fits.open(imgname)
	data = img[0].data
	img.close()
	weight_img = np.ones(data.shape)
	weight_img[np.isnan(data)] = 0
	weight_img[2048,:] = 0	 #Better way to flag this out 
	weight_hdu = fits.PrimaryHDU(weight_img)
	weight_hdu.writeto(imgname+'.weight',overwrite=True)


def run_swarp(imgname,weightimgname,pixscale,imgpisxize,ra,dec):
	#Run swarp on the proc image file
	try:
		#'swarp wasp_ZTF21abcwyvi24.proc.autoastrom.scampastrom.fits,wasp_ZTF21abcwyvi25.proc.autoastrom.scampastrom.fits,wasp_ZTF21abcwyvi29.proc.autoastrom.scampastrom.fits -c config.swarp -RESAMPLE Y -RESAMPLE_DIR . -SUBTRACT_BACK Y -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE wasp_ZTF21abcwyvi24.proc.autoastrom.scampastrom.fits.weight,wasp_ZTF21abcwyvi25.proc.autoastrom.scampastrom.fits.weight,wasp_ZTF21abcwyvi29.proc.autoastrom.scampastrom.fits.weight -COMBINE Y -COMBINE_TYPE MEDIAN -IMAGEOUT_NAME wasp_ZTF21abcwyvi24-25-29.autoastrom.scampastrom.fits -WEIGHTOUT_NAME wasp_ZTF21abcwyvi24-25-29.autoastrom.scampastrom.fits.weight'
		command = 'swarp ' + imgname + ' -c ' + astrom_swarp + '  -RESAMPLE Y -RESAMPLE_DIR ' + '.' + ' -SUBTRACT_BACK N -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE ' + weightimgname + ' ' + ' -COMBINE Y -COMBINE_TYPE MEDIAN -PIXEL_SCALE %.3f'%(pixscale) + ' -IMAGE_SIZE %i'%(imgpisxize) + ' -CENTER %.8f,%.8f'%(ra,dec) + ' -COPY_KEYWORDS FILTER,OBJECT ' + ' -IMAGEOUT_NAME ' + imgname.replace('.fits','') + '.resamp.fits' +' -WEIGHTOUT_NAME ' + imgname.replace('.fits','') + '.resamp.weight'
		print('Executing command : %s'%(command))
		rval = subprocess.run(command.split(),check=True,capture_output=True)
		print('Process completed')
		print(rval.stdout.decode())
		return 0

	except subprocess.CalledProcessError as err:
		print('Could not run swarp with error %s.'%(err))
		return -1

def firstPass_astrometry(procimgname,pixscale=0.466):
	#Run autoastrometry.py
	try:
		command = 'python autoastrometry.py ' + procimgname + ' -px '+ str(pixscale) + ' -inv -o ' + procimgname.replace('.fits','') + '.autoastrom.fits' #Add additional parameters based on pixelscale etc.
		print('Executing command : %s'%(command))
		rval = subprocess.run(command.split(),check=True,capture_output=True)
		print('Process completed')
		print(rval.stdout.decode())
		return 0

	except subprocess.CalledProcessError as err:
		print('Could not run astrometry with error %s.'%(err))
		return -1


def run_astrometry(proclist,pixscale=0.466,weightimg='weight.fits'):
	for procimgname in proclist:
		res = firstPass_astrometry(procimgname,pixscale)
		if res<0:
			print('Failed at first pass (autoastrometry) stage for ',procimgname)
			continue
		
		res = secondPass_astrometry(procimgname.replace('.fits','') + '.autoastrom.fits',writefiles=True,weightimg=weightimg)
		if res<0:
			print('Failed at second pass stage for ',procimgname)
			continue
		
		gen_weight_img(procimgname.replace('.fits','')+'.autoastrom.scampastrom.fits')
		img = fits.open(procimgname.replace('.fits','')+'.autoastrom.fits')
		header = img[0].header
		data = img[0].data
		img.close()

		w = WCS(header)
		cd11 = header['CD1_1']
		cd21 = header['CD2_1']
		cd12 = header['CD1_2']
		cd22 = header['CD2_2']
		nxpix = header['NAXIS1']
		nypix = header['NAXIS2']

		image_x_cen = nxpix/2
		image_y_cen = nypix/2

		[raCent, decCent] = w.all_pix2world(image_x_cen,image_y_cen,1)

		xscale = np.sqrt(cd11**2+cd21**2)
		yscale = np.sqrt(cd12**2 + cd22**2)

		if not os.path.exists('resamp'):
			os.mkdir('resamp')
			print('Had to make directory resamp')
		run_swarp(procimgname.replace('.fits','')+'.autoastrom.scampastrom.fits',procimgname.replace('.fits','')+'.autoastrom.scampastrom.fits'+'.weight',xscale*3600,max(nxpix,nypix)*np.sqrt(2),raCent,decCent)


if __name__ =='__main__':
	import argparse
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--d",type=str,default='data',help='Path to the directory where raw images are stored')
	parser.add_argument("--outdir",type=str,default='proc',help='Name of output subdirectory for processed images')
	parser.add_argument("--writefiles",action="store_false",help='Write the new scamped astrometric files?')
	parser.add_argument("--redo",action="store_true",help='Redo astrometry')
	parser.add_argument("--single",action="store_true",help='Single image')
	parser.add_argument("--imagename",type=str,help='Single image')
	parser.add_argument("--pixscale",default=0.466,help='Pixelscale')
	parser.add_argument("--weightimg",type=str,help='weight image')
	
	
	args = parser.parse_args()

	if args.single:
		run_astrometry([args.imagename],args.pixscale,args.weightimg)
		exit(0)

	proclist = glob('%s/%s/*.fits'%(args.d,args.outdir))
	allproclist = proclist.copy()
	for procfile in proclist:
		if os.path.exists('%s.autoastrom.scampastrom.resamp.fits'%(procfile.replace('.fits',''))) and not args.redo:
			proclist.remove(procfile)


	print('Found %s total files, %s of which have no astrometry.'%(len(allproclist),len(proclist)))
	run_astrometry(proclist)


