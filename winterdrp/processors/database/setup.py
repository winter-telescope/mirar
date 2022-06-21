import psycopg
from astropy.io import fits
import os
from glob import glob
from astropy.time import Time
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from winterdrp.processors.database.paths import schema_dir
from psycopg.errors import Error
import logging

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

def check_if_table_exists(db_name, db_table, db_user, password):  # db_host?
    with psycopg.connect(f"user={db_user} password={password}") as conn:
        conn.autocommit = True
        command = '''SELECT datname FROM pg_database;'''
        data = conn.execute(command).fetchall()

    existing_db_names = [x[0] for x in data]
    logger.debug(f"Found the following databases: {existing_db_names}")

    db_exist_bool = db_name in existing_db_names
    logger.info(f"Database table '{db_table}' {['does not exist', 'already exists'][db_exist_bool]}")

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
    header, db_user, password, db_name, db_table
):
    with psycopg.connect(f"dbname={db_name} user={db_user} password={password}") as conn:
        conn.autocommit = True

        sql_query = f"""
        SELECT Col.Column_Name from
            INFORMATION_SCHEMA.TABLE_CONSTRAINTS Tab,
            INFORMATION_SCHEMA.CONSTRAINT_COLUMN_USAGE Col
        WHERE
            Col.Constraint_Name = Tab.Constraint_Name
            AND Col.Table_Name = Tab.Table_Name
            AND Constraint_Type = 'PRIMARY KEY'
            AND Col.Table_Name = '{db_table}'
        """

#         sql_query = f"""
#
# SELECT COLUMN_NAME
# FROM INFORMATION_SCHEMA.COLUMNS
# WHERE TABLE_NAME = {db_table}
# EXCEPT
# SELECT COLUMN_NAME
# FROM INFORMATION_SCHEMA.TABLE_CONSTRAINTS [tc]
# JOIN INFORMATION_SCHEMA.KEY_COLUMN_USAGE [ku] ON tc.CONSTRAINT_NAME = ku.CONSTRAINT_NAME AND ku.table_name = {db_table} AND tc.CONSTRAINT_TYPE = 'PRIMARY KEY'
# """

        # primary_key = conn.execute(sql_query).fetchone()
        # print(primary_key)

        # sql_query = f"""
        #
        # """

        # with conn.execute(f"SELECT * FROM {db_table} LIMIT 0") as cursor:
        with conn.execute(sql_query) as cursor:

            primary_key = [x[0] for x in cursor.fetchall()]

            colnames = [
                desc[0] for desc in conn.execute(f"SELECT * FROM {db_table} LIMIT 1").description
                if desc[0] not in primary_key
            ]

            logger.debug(colnames)

            txt = f'INSERT INTO {db_table} ({colnames}) VALUES ('

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
            export_to_db(header, db_name=db_name, db_user=user, password=pwd)
        except Error as e:
            err = f"Error processing image {path} :'{e}'"
            logger.error(err)
        # except FileExistsError:
        #     pass
                # except KeyboardInterrupt:
                #     pass

