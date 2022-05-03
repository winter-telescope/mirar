import psycopg2
from astropy.coordinates import SkyCoord
import astropy.units as u
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from astropy.io import fits, ascii
import os
from glob import glob
from astropy.time import Time


def ingest_raw(user,password,rawimlist,reingest=False):
	conn = psycopg2.connect(database='commissioning', user=user, password=password)#, host='127.0.0.1', port='5432')
	conn.autocommit = True
	cursor = conn.cursor()
	for l in rawimlist:

		#Check if entry exists before
		command = '''SELECT raw.rawid FROM raw WHERE raw.filename = '%s' '''%(l)
		cursor.execute(command)
		match_id = cursor.fetchall()
		if len(match_id)>0:
			if not reingest:
				continue
			if reingest:
				raw_matchid = match_id[0][0]
				command = '''DELETE FROM raw WHERE raw.rawid = %i'''%(raw_matchid)
				cursor.execute(command)
				print('Deleted 1 row')

		img = fits.open(l)
		header = img[0].header
		img.close()

		filename = l

		if 'UTCISO' in header.keys():
			obs_utc = header['UTCISO']
			obsdate = int(Time(obs_utc).isot.split('T')[0].replace('-',''))

		
		filesize = os.stat(filename).st_size
		
		if 'OBJECT' in header.keys():
			object_name = header['OBJECT']
		else:
			object_name = ''

		if 'OBSTYPE' in header.keys():
			obstype = header['OBSTYPE']
		else:
			obstype = ''

		if 'RA' in header.keys():
			obj_ra = header['RA']
			obj_dec = header['DEC']
			obj_crds = SkyCoord(ra=obj_ra,dec=obj_dec,unit=(u.deg,u.deg))
			obj_ra, obj_dec = obj_crds.ra.deg, obj_crds.dec.deg
		else:
			obj_ra, obj_dec = -99,-99

		if 'TELRA' in header.keys():
			tel_ra = header['TELRA']
			tel_dec = header['TELDEC']
			tel_crds = SkyCoord(ra=tel_ra,dec=tel_dec,unit=(u.deg,u.deg))
			tel_ra, tel_dec = tel_crds.ra.deg, tel_crds.dec.deg
		else:
			tel_ra, tel_dec = -99,-99

		if 'HA' in header.keys():
			ha = header['HA']
		else:
			ha = -99

		if 'AIRMASS' in header.keys():
			airmass = header['AIRMASS']
		else:
			airmass = -99

		if 'ALTITUDE' in header.keys():
			alt = header['ALTITUDE']
			az = header['AZIMUTH']
		else:
			alt,az = -99, -99

		if 'MOONANGLE' in header.keys():
			moon_angle = header['MOONANGLE']
		else:
			moon_angle = -99

		if 'FIELDID' in header.keys():
			if header['FIELDID'] == '':
				field_id = -99
			else:
				field_id = header['FIELDID']
		else:
			field_id = -99

		field_id = -99
		#print(field_id)

		if 'FOCPOS' in header.keys():
			focuser_position = header['FOCPOS']
		else:
			focuser_position = -99

		if 'ROTFIELD' in header.keys():
			rot_field_angle = header['ROTFIELD']
		else:
			rot_field_angle = -999

		if 'ROTMECH' in header.keys():
			rot_position = header['ROTMECH']
		else:
			rot_position = -999

		if 'DOME_AZ' in header.keys():
			dome_az = header['DOME_AZ']
		else:
			dome_az = -999

		if 'DOME_SHUTTER' in header.keys():
			dome_shutter = header['DOME_SHUTTER']
		else:
			dome_shutter = ''

		if 'FILTERID' in header.keys():
			filt = header['FILTERID'][:10]
			if 'other' in filt:
				filt = 'r'
		else:
			filt = ''


		if 'SHUTOPEN' in header.keys():
			exp_start = header['SHUTOPEN']
			try:
				exp_start_mjd = Time(exp_start).mjd
			except ValueError:
				exp_start_mjd = -99
		else:
			exp_start_mjd = -99



		if 'EXPTIME' in header.keys():
			exptime = header['EXPTIME']
		else:
			exptime = -99


		if 'CCD_TEMP' in header.keys():
			ccd_temp = header['CCD_TEMP']
		else:
			ccd_temp = -99


		if 'GAIN' in header.keys():
			ccd_gain = header['GAIN']
		else:
			ccd_gain = -99


		if 'FREQ' in header.keys():
			ccd_freq = header['FREQ']
		else:
			ccd_freq = -99
		
		progid = 1
		camera = 'SUMMER'
		#telescop = 'WINTER'
		print(filename,obsdate,filesize,object_name,obj_ra,obj_dec,tel_ra, tel_dec,ha, airmass, alt, az, moon_angle, field_id, focuser_position, rot_field_angle, rot_position, dome_az, dome_shutter, filt, exp_start_mjd, exptime, ccd_temp, ccd_gain, ccd_freq, camera, obstype, progid)
		command = '''INSERT INTO raw (FILENAME, OBSDATE, FILESIZE, OBJECT, OBJ_RA, OBJ_DEC, TEL_RA, TEL_DEC, HA, AIRMASS, ALT, AZ, MOON_ANGLE, FIELD_ID, FOCUSER_POSITION, ROT_FIELD_ANGLE,  ROT_POSITION, DOME_AZ, DOME_SHUTTER, FILTER, EXP_START_MJD, EXPTIME, CCD_TEMP, CCD_GAIN, CCD_READOUT_FREQ, CAMERA, OBSTYPE, PROGID) VALUES ('%s',%i,%.2f, '%s',%.5f, %.5f,%.5f, %.5f,%.2f, %.3f, %.2f, %.2f, %.3f, %i, %.2f, %.3f, %.3f, %.2f, '%s','%s', %.4f, %.2f,%.1f, %.2f, %.2f, '%s', '%s',%i);'''%(filename,obsdate,filesize,object_name,obj_ra,obj_dec,tel_ra, tel_dec,ha, airmass, alt, az, moon_angle, field_id, focuser_position, rot_field_angle, rot_position, dome_az, dome_shutter, filt, exp_start_mjd, exptime, ccd_temp, ccd_gain, ccd_freq, camera, obstype, progid)
		cursor.execute(command)
		#'''CREATE TABLE raw ( rawid SERIAL PRIMARY KEY, FILENAME VARCHAR(255) NOT NULL, OBSDATE INT, FILESIZE REAL, OBJECT VARCHAR(30), OBJ_RA FLOAT, OBJ_DEC FLOAT, TEL_RA FLOAT, TEL_DEC FLOAT, HA FLOAT, AIRMASS FLOAT, ALT FLOAT, AZ FLOAT, MOON_ANGLE FLOAT, FIELD_ID VARCHAR(30), FOCUSER_POSITION FLOAT, ROT_FIELD_ANGLE FLOAT, ROT_POSITION FLOAT, DOME_AZ FLOAT, DOME_SHUTTER VARCHAR(5), FILTER VARCHAR(5), EXP_START_MJD FLOAT, EXPTIME FLOAT, CCD_TEMP FLOAT, CCD_GAIN FLOAT, CCD_READOUT_FREQ FLOAT  );'''
		print('Ingested %s'%(l))
	return 0	


def ingest_proc(user,password,procimlist,reingest=False):
	conn = psycopg2.connect(database='commissioning', user=user, password=password)#, host='127.0.0.1', port='5432')
	conn.autocommit = True
	cursor = conn.cursor()

	for l in procimlist:
		img = fits.open(l)
		header = img[0].header
		img.close()
		proc_filename = l

		#Check if entry exists previously
		command = '''SELECT proc.procid FROM proc WHERE proc.proc_imagename = '%s' '''%(proc_filename)
		cursor.execute(command)
		match_id = cursor.fetchall()
		if len(match_id)>0:
			if not reingest:
				continue
			if reingest:
				proc_matchid = match_id[0][0]
				command = '''DELETE FROM proc WHERE proc.procid = %i'''%(proc_matchid)
				cursor.execute(command)
				print('Deleted 1 row')

		raw_filename = proc_filename.replace('.proc.autoastrom.scampastrom.resamp.fits','.fits').replace('/proc','')
		command = '''SELECT raw.rawid FROM raw WHERE raw.filename = '%s' '''%(raw_filename)
		cursor.execute(command)
		rawid = cursor.fetchall()[0][0]

		if 'UTCISO' in header.keys():
			obs_utc = header['UTCISO']
			obsdate = int(Time(obs_utc).isot.split('T')[0].replace('-',''))

		if 'FILTERID' in header.keys():
			filt = header['FILTERID'][:10]
		else:
			filt = ''

		if 'EXPTIME' in header.keys():
			exptime = header['EXPTIME']
		else:
			exptime = -99

		if 'FIELDID' in header.keys():
			field_id = header['FIELDID']
			if field_id == '':
				field_id = -99
		else:
			field_id = -99
		
		field_id = -99

		if 'FWHM_MED' in header.keys():
			fwhm = header['FWHM_MED']
		else:
			fwhm = -99

		if 'FLAT' in header.keys():
			flatfilename = header['FLAT']
		else:
			flatfilename = ''

		if 'BIAS' in header.keys():
			biasfilename = header['BIAS']
		else:
			biasfilename = ''

		if 'DARK' in header.keys():
			darkfilename = header['DARK']
		else:
			darkfilename = ''


		if 'LIMMAG' in header.keys():
			maglim = header['LIMMAG']
		else:
			maglim = -99

		#obj = header['OBJECT']
		camera = 'SUMMER'
		
		cd1_1 = header['CD1_1']
		cd1_2 = header['CD1_2']
		cd2_1 = header['CD2_1']
		cd2_2 = header['CD2_2']
		crval1 = header['CRVAL1']
		crval2 = header['CRVAL2']
		crpix1 = header['CRPIX1']
		crpix2 = header['CRPIX2']
		astr_dpa = header['ASTR_DPA']
		astr_off = header['ASTR_OFF']
		if 'FLXSCALE' in header.keys():
			flxscale = header['FLXSCALE']
		else:
			flxscale = -99
		
		pv1_0, pv1_1, pv1_4, pv1_7, pv2_0, pv2_1, pv2_2, pv2_6, pv2_10 = -99, -99, -99, -99, -99, -99, -99, -99, -99
		if 'PV1_0' in header.keys():
			pv1_0 = header['PV1_0']
		if 'PV1_1' in header.keys():
			pv1_1 = header['PV1_1']
		if 'PV1_4' in header.keys():
			pv1_4 = header['PV1_4']
		if 'PV1_7' in header.keys():
			pv1_7 = header['PV1_7']
		if 'PV2_0' in header.keys():
			pv2_0 = header['PV2_0']
		if 'PV2_1' in header.keys():
			pv2_1 = header['PV2_1']
		if 'PV2_2' in header.keys():
			pv2_2 = header['PV2_2']
		if 'PV2_6' in header.keys():
			pv2_6 = header['PV2_6']
		if 'PV2_10' in header.keys():
			pv2_10 = header['PV2_10']

		if 'ZP_AUTO' in header.keys():
			zp = header['ZP_AUTO']
		else:
			zp = -99

		if 'ZP_AUTO_NSTARS' in header.keys():
			zp_nstars = header['ZP_AUTO_NSTARS']
		else:
			zp_nstars = -99

		if 'ZP_AUTO_STD' in header.keys():
			zp_std = header['ZP_AUTO_STD']
		else:
			zp_std= -99
		
		if 'SHUTOPEN' in header.keys():
			exp_start = header['SHUTOPEN']
			exp_start_mjd = Time(exp_start).mjd
		else:
			exp_start_mjd = -99

		#resamp_img = l.replace('.proc.autoastrom.scampastrom.fits','.fits')
		#raw_l = l1.replace('proc/','')
		#zp = header['ZP_AUTO']
		
		#img = fits.open(raw_l)
		#header = img[0].header
		#img.close()
		#obsdate = header['UTCSTART']
		#object_name = header['OBJECT']
		#exposure_length = header['EXPTIME']

		command = '''INSERT INTO proc (rawid, RAW_IMAGENAME, PROC_IMAGENAME, OBSDATE, FILTER,  EXPTIME, FIELD_ID, CD1_1, CD1_2, CD2_1, CD2_2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, ZP, FWHM, MAGLIM, FLATFILENAME, DARKFILENAME, BIASFILENAME, ASTR_DPA, ASTR_OFF, FLXSCALE, PV1_0, PV1_1, PV1_4, PV1_7, PV2_0, PV2_1, PV2_2, PV2_6, PV2_10, ZP_NSTARS, ZP_STD, EXP_MJD) VALUES (%i, '%s', '%s', %i,'%s',%.3f,%i,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%.2f,%.2f,'%s','%s','%s',%.2f,%.2f,%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%i,%.4f,%.4f);'''%(rawid,raw_filename,proc_filename,obsdate,filt,exptime,field_id,cd1_1,cd1_2,cd2_1,cd2_2,crval1,crval2,crpix1,crpix2,zp,fwhm,maglim,flatfilename,darkfilename,biasfilename,astr_dpa,astr_off,flxscale,pv1_0,pv1_1,pv1_4,pv1_7,pv2_0,pv2_1,pv2_2,pv2_6,pv2_10,zp_nstars,zp_std,exp_start_mjd)
		cursor.execute(command)
		print('Ingested %s'%(l))
	return 0


def create_db(user,password,dbname='commissioning'):
	conn = psycopg2.connect(database='postgres', user=user, password=password)#, host='127.0.0.1', port='5432')
	cursor = conn.cursor()
	conn.autocommit = True
	sql = '''CREATE database %s'''%(dbname)
	cursor.execute(sql)
	print('Created db %s'%(dbname))
	
	
def create_tables(user,password,dbname='commissioning'):
	conn = psycopg2.connect(database=dbname, user=user, password=password)#, host='127.0.0.1', port='5432')
	cursor = conn.cursor()
	conn.autocommit = True
	
	command = '''CREATE TABLE raw ( rawid SERIAL PRIMARY KEY, FILENAME VARCHAR(255) NOT NULL, OBSDATE INT, 
	FILESIZE REAL, OBJECT VARCHAR(30), OBJ_RA FLOAT, OBJ_DEC FLOAT, TEL_RA FLOAT, TEL_DEC FLOAT, HA FLOAT, 
	AIRMASS FLOAT, ALT FLOAT, AZ FLOAT, MOON_ANGLE FLOAT, FIELD_ID VARCHAR(30), FOCUSER_POSITION FLOAT, 
	ROT_FIELD_ANGLE FLOAT, ROT_POSITION FLOAT, DOME_AZ FLOAT, DOME_SHUTTER VARCHAR(5), FILTER VARCHAR(5), 
	EXP_START_MJD FLOAT, EXPTIME FLOAT, CCD_TEMP FLOAT, CCD_GAIN FLOAT, CCD_READOUT_FREQ FLOAT, CAMERA VARCHAR(15), 
	OBSTYPE VARCHAR(15)  ); '''
	cursor.execute(command)

	print('Created table raw in ', dbname)

	command = '''CREATE TABLE proc ( procid SERIAL PRIMARY KEY, rawid SERIAL, RAW_IMAGENAME VARCHAR(255), 
	PROC_IMAGENAME VARCHAR(255), OBSDATE INT, FILTER VARCHAR(10),  EXPTIME INT, FIELD_ID VARCHAR(30), CD1_1 REAL, 
	CD1_2 REAL, CD2_1 REAL, CD2_2 REAL, CRVAL1 REAL, CRVAL2 REAL, CRPIX1 REAL, CRPIX2 REAL, ZP REAL, FWHM REAL, 
	MAGLIM REAL, FLATFILENAME VARCHAR(255), DARKFILENAME VARCHAR(255), BIASFILENAME VARCHAR(255), ASTR_DPA FLOAT, 
	ASTR_OFF FLOAT, FLXSCALE FLOAT, PV1_0 FLOAT, PV1_1 FLOAT, PV1_4 FLOAT, PV1_7 FLOAT, PV2_0 FLOAT, PV2_1 FLOAT, 
	PV2_2 FLOAT, PV2_6 FLOAT, PV2_10 FLOAT, ZP_NSTARS INT, ZP_STD FLOAT, EXP_MJD FLOAT);''' #add which flats were used
	# and which darks were used
	cursor.execute(command)
	print('Created table proc in ', dbname)


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--d",type=str,default='data',help='Path to the directory where raw images are stored')
	parser.add_argument("--procdir",type=str,default='proc',help='Name of output subdirectory for processed images')
	parser.add_argument("--create_db",action="store_true")
	parser.add_argument("--raw",action="store_true")
	parser.add_argument("--proc",action="store_true")
	parser.add_argument("--reingest",action="store_true")
	parser.add_argument("--dbname",type=str,default='commissioning',help='Name of database')
	parser.add_argument("--secrets",type=str,default='/home/viraj/db_secrets.csv',help='Username and password file')


	args = parser.parse_args()
	
	secrets = ascii.read(args.secrets)
	user = secrets['user'][0]
	password = secrets['pwd'][0]

	if args.create_db:
		create_db(user,password,args.dbname)

		create_tables(user,password,args.dbname)

	if args.raw:
		rawimlist = glob('%s/*.fits'%(args.d))
		print('Found %s raw files'%(len(rawimlist)))
		ingest_raw(user, password, rawimlist,reingest=args.reingest)

	if args.proc:
		procimlist = glob('%s/proc/*.resamp.fits'%(args.d))
		print('Found %s processed (scamped) files'%(len(procimlist)))
		ingest_proc(user, password, procimlist,reingest=args.reingest)
