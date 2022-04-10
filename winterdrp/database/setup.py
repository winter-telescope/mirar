import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from astropy.io import fits, ascii
import os
from glob import glob
from astropy.time import Time
import numpy as np
import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u


def create_db(db_name, db_user, password):
    conn = psycopg2.connect(database='postgres', user=db_user, password=password)  # , host='127.0.0.1', port='5432')
    cursor = conn.cursor()
    conn.autocommit = True
    sql = f'''CREATE database {db_name}'''
    cursor.execute(sql)
    print(f'Created db {db_name}')


def create_table(schema_path, db_name, db_user, password, db_host='localhost'):
    conn = psycopg2.connect(database=db_name,  user=db_user, password=password)
    conn.autocommit = True
    cursor = conn.cursor()

    cursor.execute(open(schema_path, "r").read())


def check_if_db_exists(db_name, db_user, password, db_host='localhost'):
    conn = psycopg2.connect(database='postgres',  user=db_user, password=password)
    conn.autocommit = True
    cursor = conn.cursor()

    command = '''SELECT datname FROM pg_database;'''
    cursor.execute(command)
    data = cursor.fetchall()

    conn.close()
    existing_db_names = [x[0] for x in data]
    # print(existing_db_names)
    return db_name in existing_db_names


def insert_exposures(header, db_user, password, db_name='commissioning', db_host='localhost'):
    conn = psycopg2.connect(database=db_name, user=db_user, password=password)
    conn.autocommit = True
    cursor = conn.cursor()

    cursor.execute("SELECT * FROM exposures LIMIT 0")
    colnames = [desc[0] for desc in cursor.description]
    print(colnames)
    colnames.remove('expid')
    #header['EXPID'] = 0
    header['OBSDATE'] = int(header['UTC'].split('_')[0])
    obstime = Time(header['UTCISO'],format='iso')
    t0 = Time('2018-01-01',format='iso')
    header['OBSHISTID'] = 0
    header['NIGHT'] = np.floor(((obstime - t0)).jd
                                   ).astype(np.int)
    header['PROGID'] = 0
    header['VISITTIME'] = header['AEXPTIME']
    header['VISITEXPTIME'] = header['AEXPTIME']
    header['EXPMJD'] = header['OBSMJD']
    header['SUBPROGRAM'] = 'high_cadence'
    header['FILTER'] = header['FILTERID']

    print('Time',header['shutopen'])
    try:
        header['shutopen'] = Time(header['SHUTOPEN'],format='iso').jd
    except ValueError:
        header['shutopen'] = -99

    try:
        header['shutclsd'] = Time(header['SHUTCLSD'], format='iso').jd
    except:
        header['shutclsd'] = -99
    header['PROCESS_FLAG'] = 0
    sunmoon_keywords = ['Moonra', 'Moondec', 'Moonillf', 'Moonphas', 'Moonalt', 'Sunaz', 'Sunalt']
    for key in sunmoon_keywords:
        if header[key] == '':
            header[key] = 0


    itid_dict = {'SCIENCE':1,
               'BIAS':2, 'FLAT':2, 'DARK':2, 'FOCUS' :3
               ,'POINTING':4, 'OTHER':5}

    header['ITID'] = itid_dict[header['OBSTYPE']]

    if header['FIELDID'] == 'radec':
        header['FIELDID'] = 0

    if header['ITID'] != 1:
        header['FIELDID'] = -99


    crds = SkyCoord(ra=header['RA'],dec=header['DEC'],unit=(u.deg,u.deg))
    header['RA'] = crds.ra.deg
    header['DEC'] = crds.dec.deg

    txt = f'INSERT INTO exposures ({colnames}) VALUES ('
    txt = txt.replace('[','').replace(']','').replace("'",'')

    for c in colnames:
        txt = txt + "'" + str(header[c]) + "'"+ ', '
    txt = txt+');'
    txt = txt.replace(', )',')')
    print(txt)
    command = txt
    cursor.execute(command)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', type=str, default='commissioning_1', help='Database name')
    parser.add_argument('--secrets', type=str, default='/Users/viraj/winter_drp/db_secrets.csv')
    parser.add_argument('--h', type=str, default='localhost')
    parser.add_argument('--d',type=str,default='/Users/viraj/winter_data/summer/20220402/raw/')

    args = parser.parse_args()

    db_name = args.db
    secrets_file = args.secrets
    host = args.h
    dirname = args.d

    secrets = ascii.read(secrets_file)
    user = secrets['user'][0]
    pwd = secrets['pwd'][0]

    if not check_if_db_exists(db_name, user, pwd, host):
        create_db(db_name, user, pwd)

    table_names = ['fields', 'filters', 'nights', 'itid', 'exposures', 'programs']#, 'detfiles']
    for table in table_names:
        create_table('schema/' + table+'.sql', db_name, user, pwd)


    ls = glob(os.path.join(dirname, 'SUMMER_*.fits'))
    print(len(ls))
    #header = fits.getheader('/Users/viraj/winter_data/summer/20220402/raw/SUMMER_20220402_203525_Camera0.fits')
    for l in ls:
        header = fits.getheader(l)
        print(l)
        insert_exposures(header, db_name=db_name, db_user=user, password=pwd, db_host=host)

