from glob import glob
from astropy.io import fits
import numpy as np
import argparse

dates = [20210802,20210804,20210805,20210806,20210809,20210810]

for date in dates:
	ls = glob('/scr2/viraj/winter_data/commissioning/raw/%s/*.fits'%(date))
	ls = np.sort(ls)
	with open('/scr2/viraj/winter_data/commissioning/raw/%s/%s_log.csv'%(date,date),'w') as f:
	    f.write('Name,RA,Dec,exptime,median_counts\n')
	    for l in ls:
	        img = fits.open(l)
	        header = img[0].header
	        data = img[0].data
	        img.close()
	        f.write('%s,%s,%s,%s,%.2f\n'%(l,header['RA'],header['DEC'],header['AEXPTIME'],np.median(data)))
