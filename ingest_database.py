import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from astropy.io import fits
import os
from glob import glob


def ingest_raw(rawimlist):
	conn = psycopg2.connect(database='testdb3', user='viraj', password='', host='127.0.0.1', port='5432')
	conn.autocommit = True
	cursor = conn.cursor()
	for l in rawimlist:
		img = fits.open(l)
		header = img[0].header
		img.close()
		image_name = l
		obsdate = header['UTCSTART']
		filename = l
		filesize = os.stat(filename).st_size
		filt = header['FILTER'][0]
		exposure_length = header['EXPTIME']
		object_name = header['OBJECTA']
		object_type = header['OBJECT']
		telescop = 'WINTER'
		command = '''INSERT INTO raw (image_name, OBSDATE, FILENAME, FILESIZE, FILTER, EXPOSURE_LENGTH, OBJECT, TELESCOP) VALUES ('%s','%s','%s',%i,'%s',%i,'%s','%s');'''%(l,obsdate,l,filesize,filt,exposure_length,object_name,telescop)
		cursor.execute(command)

	return 0	


def ingest_proc(procimlist):
	conn = psycopg2.connect(database='testdb3', user='viraj', password='', host='127.0.0.1', port='5432')
	conn.autocommit = True
	cursor = conn.cursor()

	for l in procimlist:
		l = ls[0]
		img = fits.open(l)
		header = img[0].header
		img.close()
		image_name = l
		filt = header['FILTER'][0]
		obj = 'ZTF21abcwyvi'#header['OBJECT']
		telescop = 'WINTER'
		cd1_1 = header['CD1_1']
		cd1_2 = header['CD1_2']
		cd2_1 = header['CD2_1']
		cd2_2 = header['CD2_2']
		crval1 = header['CRVAL1']
		crval2 = header['CRVAL2']
		crpix1 = header['CRPIX1']
		crpix2 = header['CRPIX2']
		zp = header['HIERARCH ZP_6.0']
		l1 = l.replace('.proc.autoastrom.scampastrom.resamp.fits','.fits')
		raw_l = l1.replace('proc/','')

		command = '''SELECT raw.rawid FROM raw WHERE raw.image_name = '%s' '''%(raw_l)
		cursor.execute(command)
		img = fits.open(raw_l)
		header = img[0].header
		img.close()
		obsdate = header['UTCSTART']
		object_name = header['OBJECT']
		exposure_length = header['EXPTIME']

		rawid = cursor.fetchall()[0][0]


def create_db(dbname='commissioning'):
	conn = psycopg2.connect(database='postgres', user='viraj', password='', host='127.0.0.1', port='5432')
	cursor = conn.cursor()
	conn.autocommit = True
	sql = '''CREATE database %s'''%(dbname)
	cursor.execute(sql)
	print('Created db %s'%(dbname))
	
	
def create_tables(dbname='commissioning'):
	conn = psycopg2.connect(database=dbname, user='viraj', password='', host='127.0.0.1', port='5432')
	cursor = conn.cursor()
	conn.autocommit = True
	
	command = '''CREATE TABLE raw ( rawid SERIAL PRIMARY KEY, FILENAME VARCHAR(255) NOT NULL, OBSDATE INT, FILESIZE REAL, OBJECT VARCHAR(30), OBJ_RA FLOAT, OBJ_DEC FLOAT, TEL_RA FLOAT, TEL_DEC FLOAT, HA FLOAT, AIRMASS FLOAT, ALT FLOAT, AZ FLOAT, MOON_ANGLE FLOAT, FIELD_ID VARCHAR(30), FOCUSER_POSITION FLOAT, ROT_FIELD_ANGLE FLOAT, ROT_POSITION FLOAT, DOME_AZ FLOAT, DOME_SHUTTER VARCHAR(5), FILTER VARCHAR(5), EXP_START_MJD FLOAT, EXPTIME FLOAT, CCD_TEMP FLOAT, CCD_GAIN FLOAT, CCD_READOUT_FREQ FLOAT  );'''
	cursor.execute(command)

	print('Created table raw in ', dbname)

	command = '''CREATE TABLE proc ( procid SERIAL PRIMARY KEY, rawid SERIAL, RAW_IMAGENAME VARCHAR(255), PROC_IMAGENAME VARCHAR(255), OBSDATE INT, FILTER VARCHAR(10),  EXPTIME INT, FIELD_ID VARCHAR(30), CD1_1 REAL, CD1_2 REAL, CD2_1 REAL, CD2_2 REAL, CRVAL1 REAL, CRVAL2 REAL, CRPIX1 REAL, CRPIX2 REAL, ZP REAL, FWHM REAL, MAGLIM REAL, FLATFILENAME VARCHAR(255), DARKFILENAME VARCHAR(255), BIASFILENAME VARCHAR(255));''' #add which flats were used and which darks were used
	cursor.execute(command)
	print('Created table proc in ', dbname)


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--d",type=str,default='data',help='Path to the directory where raw images are stored')
	parser.add_argument("--procdir",type=str,default='proc',help='Name of output subdirectory for processed images')
	parser.add_argument("--create_db",action="store_true")
	parser.add_argument("--dbname",type=str,default='commissioning',help='Name of database')

	args = parser.parse_args()
	if args.create_db:
		create_db(args.dbname)

		create_tables(args.dbname)