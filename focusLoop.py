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
import sys
sys.path.append('.')
from genstats import get_img_fwhm
from scipy.optimize import curve_fit

def parabola(x,x0,A,B):
	return A + B*(x-x0)**2


def fit_parabola(focus_vals,fwhms,stds,plot=True):
	p0 = [np.mean(focus_vals),np.min(fwhms),np.std(fwhms)]
	popt,pcov = curve_fit(parabola,xdata=focus_vals,ydata=fwhms,p0=p0,sigma=stds)
	if plot:
		plt.figure()
		plt.errorbar(focus_vals,fwhms,yerr=stds,fmt='.',c='red')
		plotfoc = np.linspace(np.min(focus_vals),np.max(focus_vals),20)
		print(popt)
		plt.plot(plotfoc,parabola(plotfoc,popt[0],popt[1],popt[2]))
		plt.title('Best FWHM : %.1f arcsec'%(np.min(fwhms)))
		plt.savefig('focusloop.pdf',bbox_inches='tight')
	return popt


def analyse_imgs_focus(imglist,pixscale=0.466,exclude=False):
	med_fwhms = []
	std_fwhms = []
	focus_vals = []
	for imgname in imglist:
		print(imgname)
		img = fits.open(imgname)
		header = img[0].header
		img.close()
		print(header['focuser_position']%50)
		if header['focuser_position']%50 < 5 or header['focuser_position']%50 >45:
			mean, med, std = get_img_fwhm(imgname,pixscale=pixscale,exclude=exclude)
			med_fwhms.append(med*pixscale)
			std_fwhms.append(std*pixscale)
			focus_vals.append(header['focuser_position'])
			print(med)

		#focus_vals.append()

	med_fwhms = np.array(med_fwhms)
	std_fwhms = np.array(std_fwhms)
	focus_vals = np.array(focus_vals)

	print(med_fwhms,std_fwhms,focus_vals)
	best_pars = fit_parabola(focus_vals,med_fwhms,std_fwhms)

	return best_pars[0]


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--d",type=str,help='Dirname')

	parser.add_argument("--pixscale",type=float,default=0.466,help='Pixelscale')
	parser.add_argument("--plot",action='store_true',help='Plot?')	
	parser.add_argument("--exclude",action='store_true',help='Exclude central values')	
	args = parser.parse_args()

	filterlist = ['u','r']
	filt_numlist = {'1':'u','3':'r'}
	
	focuslist = {}
	focuslist['u'] = []
	focuslist['r'] = []

	imglist = glob('%s/*.fits'%(args.d))
	for imgname in imglist:
		img = fits.open(imgname)
		header = img[0].header
		img.close()

		if header['OBSTYPE'] == 'focus':
			imfilt = str(header['FILTER'])
			
			focuslist[imfilt].append(imgname)


	for filt in filterlist:
		print('Found %s focus files for filter %s'%(len(focuslist[filt]),filt))
		if len(focuslist[filt])==0:
			print('Skipping')
			continue
		bestfoc = analyse_imgs_focus(focuslist[filt],exclude=args.exclude)
		print(bestfoc)
