import psycopg
from astropy.io import fits, ascii
import os
from glob import glob
from astropy.time import Time
import numpy as np
import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u
from winterdrp.database.schema import schema_dir
from psycopg.errors import Error
import logging
import traceback

logger = logging.getLogger(__name__)


def create_db(db_name, db_user, password):
    with psycopg.connect(f"user={db_user} password={password}") as conn:
        conn.autocommit = True
        sql = f'''CREATE database {db_name}'''
        conn.execute(sql)
        logger.info(f'Created db {db_name}')


def create_table(schema_path, db_name, db_user, password, db_host='localhost'):  # db_host?
    with psycopg.connect(f"dbname={db_name} user={db_user} password={password}") as conn:
        conn.autocommit = True
        with open(schema_path, "r") as f:
            conn.execute(f.read())


def check_if_db_exists(db_name, db_user, password):  # db_host?
    with psycopg.connect(f"user={db_user} password={password}") as conn:
        conn.autocommit = True
        command = '''SELECT datname FROM pg_database;'''
        data = conn.execute(command).fetchall()

    existing_db_names = [x[0] for x in data]
    logger.debug(f"Found the following databases: {existing_db_names}")

    db_exist_bool = db_name in existing_db_names
    logger.info(f"Database '{db_name}' {['does not exist', 'already exists'][db_exist_bool]}")

    return db_exist_bool


def insert_exposures(header2, db_user, password, db_name='commissioning'):

    header = dict(header2)

    header['OBSDATE'] = int(header['UTC'].split('_')[0])
    obstime = Time(header['UTCISO'], format='iso')
    t0 = Time('2018-01-01', format='iso')
    header['OBSHISTID'] = 0
    header['NIGHT'] = np.floor(((obstime - t0)).jd
                               ).astype(int)
    header['PROGID'] = 0
    header['VISITTIME'] = header['AEXPTIME']
    header['VISITEXPTIME'] = header['AEXPTIME']
    header['EXPMJD'] = header['OBSMJD']
    header['SUBPROGRAM'] = 'high_cadence'
    header['FILTER'] = header['FILTERID']

    # print('Time', header['shutopen'])
    try:
        header['SHUTOPEN'] = Time(header['SHUTOPEN'], format='iso').jd
    except (KeyError, ValueError):
        # header['SHUTOPEN'] = None
        pass

    try:
        header['SHUTCLSD'] = Time(header['SHUTCLSD'], format='iso').jd
    except ValueError:
        pass
        # header['SHUTCLSD'] = None

    header['PROCESS_FLAG'] = 0
    sunmoon_keywords = ['MOONRA', 'MOONDEC', 'MOONILLF', 'MOONPHAS', 'MOONALT', 'SUNAZ', 'SUNALT']
    for key in sunmoon_keywords:
        val = 0
        if key in header.keys():
            if header[key] not in ['']:
                val = header[key]
        header[key] = val

    itid_dict = {
        'SCIENCE': 1,
        'BIAS': 2,
        'FLAT': 2,
        'DARK': 2,
        'FOCUS': 3,
        'POINTING': 4,
        'OTHER': 5
    }

    if not header['OBSTYPE'] in itid_dict.keys():
        header['ITID'] = 5
    else:
        header['ITID'] = itid_dict[header['OBSTYPE']]

    if header['FIELDID'] == 'radec':
        header['FIELDID'] = 0

    if header['ITID'] != 1:
        header['FIELDID'] = -99

    crds = SkyCoord(ra=header['RA'], dec=header['DEC'], unit=(u.deg, u.deg))
    header['RA'] = crds.ra.deg
    header['DEC'] = crds.dec.deg

    with psycopg.connect(f"dbname={db_name} user={db_user} password={password}") as conn:
        conn.autocommit = True

        with conn.execute("SELECT * FROM exposures LIMIT 0") as cursor:
            colnames = [desc[0] for desc in cursor.description]

            logger.debug(colnames)
            colnames.remove('expid')

            txt = f'INSERT INTO exposures ({colnames}) VALUES ('

            for char in ["[", "]", "'"]:
                txt = txt.replace(char, '')

            for c in colnames:
                logger.debug(f"{c}, {header[c.upper()]}")
                txt += f"'{str(header[c.upper()])}', "

            txt = txt+');'
            txt = txt.replace(', )', ')')
            logger.debug(txt)
            command = txt
            cursor.execute(command)

def export_to_db(
    header2, fields_to_export, db_user, password, db_name='commissioning'
):
    header = dict(header2)

    # header['OBSDATE'] = int(header['UTC'].split('_')[0])
    # obstime = Time(header['UTCISO'], format='iso')
    # t0 = Time('2018-01-01', format='iso')
    # header['OBSHISTID'] = 0
    # header['NIGHT'] = np.floor(((obstime - t0)).jd
    #                            ).astype(int)
    # header['PROGID'] = 0
    # header['VISITTIME'] = header['AEXPTIME']
    # header['VISITEXPTIME'] = header['AEXPTIME']
    # header['EXPMJD'] = header['OBSMJD']
    # header['SUBPROGRAM'] = 'high_cadence'
    # header['FILTER'] = header['FILTERID']
    #
    # # print('Time', header['shutopen'])
    # try:
    #     header['SHUTOPEN'] = Time(header['SHUTOPEN'], format='iso').jd
    # except (KeyError, ValueError):
    #     # header['SHUTOPEN'] = None
    #     pass
    #
    # try:
    #     header['SHUTCLSD'] = Time(header['SHUTCLSD'], format='iso').jd
    # except ValueError:
    #     pass
    #     # header['SHUTCLSD'] = None
    #
    # header['PROCESS_FLAG'] = 0
    # sunmoon_keywords = ['MOONRA', 'MOONDEC', 'MOONILLF', 'MOONPHAS', 'MOONALT', 'SUNAZ', 'SUNALT']
    # for key in sunmoon_keywords:
    #     val = 0
    #     if key in header.keys():
    #         if header[key] not in ['']:
    #             val = header[key]
    #     header[key] = val
    #
    # itid_dict = {
    #     'SCIENCE': 1,
    #     'BIAS': 2,
    #     'FLAT': 2,
    #     'DARK': 2,
    #     'FOCUS': 3,
    #     'POINTING': 4,
    #     'OTHER': 5
    # }
    #
    # if not header['OBSTYPE'] in itid_dict.keys():
    #     header['ITID'] = 5
    # else:
    #     header['ITID'] = itid_dict[header['OBSTYPE']]
    #
    # if header['FIELDID'] == 'radec':
    #     header['FIELDID'] = 0
    #
    # if header['ITID'] != 1:
    #     header['FIELDID'] = -99
    #
    # crds = SkyCoord(ra=header['RA'], dec=header['DEC'], unit=(u.deg, u.deg))
    # header['RA'] = crds.ra.deg
    # header['DEC'] = crds.dec.deg

    with psycopg.connect(f"dbname={db_name} user={db_user} password={password}") as conn:
        conn.autocommit = True

        with conn.execute("SELECT * FROM exposures LIMIT 0") as cursor:
            colnames = [desc[0] for desc in cursor.description]

            logger.debug(colnames)
            colnames.remove('expid')

            txt = f'INSERT INTO exposures ({colnames}) VALUES ('

            for char in ["[", "]", "'"]:
                txt = txt.replace(char, '')

            for c in colnames:
                logger.debug(f"{c}, {header[c.upper()]}")
                txt += f"'{str(header[c.upper()])}', "

            txt = txt + ');'
            txt = txt.replace(', )', ')')
            logger.debug(txt)
            command = txt
            cursor.execute(command)


def set_up_db(
        data_dir: str,
        db_name: str = 'commissioning_1',
        user: str = os.path.basename(os.environ["HOME"]),
):

    # We should find a more secure method for this
    # secrets = ascii.read(secrets_file)
    # user = secrets['user'][0]
    # pwd = secrets['pwd'][0]

    # password = getpass.getpass(f"Enter password for user {user}: ")
    pwd = "nothingfornow"

    if not check_if_db_exists(db_name, user, pwd):
        create_db(db_name, user, pwd)

    table_names = ['fields', 'filters', 'nights', 'itid', 'exposures', 'programs']#, 'detfiles']
    for table in table_names:
        schema = os.path.join(schema_dir, table+'.sql')
        create_table(schema, db_name, user, pwd)

    image_list = glob(os.path.join(data_dir, 'SUMMER_*.fits'))
    for path in image_list:
        header = fits.getheader(path)
        try:
            insert_exposures(header, db_name=db_name, db_user=user, password=pwd)
            export_to_db(header, fields_to_export=[], db_name=db_name, db_user=user, password=pwd)
        except Error as e:
            err = f"Error processing image {path} :'{e}'"
            logger.error(err)
        # except FileExistsError:
        #     pass
                # except KeyboardInterrupt:
                #     pass

