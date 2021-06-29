import subprocess
import os
from datetime import datetime
from astropy.time import Time

def transfer_nightly_data():
	utc = datetime.utcnow()
	time = Time(utc)
	utc_date = time.isot.split('T')[0].replace('-','')

	if not os.path.exists('/home/winter/data/images/%s'%(utc_date)):
		return -1
	command = 'rsync -azP /home/winter/data/images/%s'%(utc_date) + ' viraj@gayatri.caltech.edu:/scr2/viraj/winter_data/commissioning/raw/'
	subprocess.call(command,shell=True)


if __name__ == '__main__':
	transfer_nightly_data()