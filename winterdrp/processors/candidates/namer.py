import logging
import os
import logging
from winterdrp.processors.base_processor import BaseDataframeProcessor
import numpy as np
import pandas as pd
from winterdrp.processors.database.postgres import execute_query, xmatch_import_db
from astropy.time import Time

logger = logging.getLogger(__name__)


class CandidateNamer(BaseDataframeProcessor):
    def __init__(self,
                 db_name: str,
                 base_name: str,
                 xmatch_radius_arcsec: float,
                 name_start: str = 'aaaaa',
                 cand_table_name: str = 'candidates',
                 db_user: str = os.environ.get('DB_USER'),
                 db_pwd: str = os.environ.get('DB_PWD'),
                 db_name_field: str = 'objectId',
                 db_order_field: str = 'candid',
                 date_field: str = 'jd',
                 *args,
                 **kwargs):
        super(CandidateNamer, self).__init__(*args, **kwargs)
        self.db_name = db_name
        self.db_user = db_user
        self.db_pwd = db_pwd
        self.cand_table_name = cand_table_name
        self.db_name_field = db_name_field
        self.db_order_field = db_order_field
        self.base_name = base_name
        self.name_start = name_start
        self.date_field = date_field
        self.xmatch_radius_arcsec = xmatch_radius_arcsec

    @staticmethod
    def increment_string(string: str):
        '''

        Parameters
        ----------
        string

        Returns
        -------
        An incremented string, eg. aaa -> aab, aaz -> aba, azz -> baa, zzz-> aaaa
        '''
        charpos = len(string) - 1
        # will iteratively try to increment characters starting from the last
        inctrue = False
        newstring = ''
        while charpos >= 0:
            cref = string[charpos]
            if inctrue:
                newstring = cref + newstring
                charpos -= 1
                continue
            creford = ord(cref)
            # increment each character, if at 'z', increment the next one
            if creford + 1 > 122:
                newstring = 'a' + newstring
                if charpos == 0:
                    newstring = 'a' + newstring
            else:
                nextchar = chr(creford + 1)
                newstring = nextchar + newstring
                inctrue = True
            charpos -= 1
            continue

        return newstring

    def get_next_name(self, candjd, lastname=None):
        candyear = Time(candjd, format='jd').datetime.year % 1000
        if lastname is None:
            query = f"""SELECT name FROM {self.cand_table_name} ORDER BY {self.db_order_field} desc LIMIT 1;"""
            res = execute_query(query, db_name=self.db_name, db_user=self.db_user, password=self.db_pwd)
            if len(res) == 0:
                name = self.base_name + str(candyear) + self.name_start
                return name
            else:
                lastname = res[0][0]
                logger.debug(res)
        lastyear = int(lastname[len(self.base_name):len(self.base_name) + 2])
        if candyear != lastyear:
            name = self.base_name + str(candyear) + self.name_start
        else:
            lastname_letters = lastname[len(self.base_name) + 2:]
            newname_letters = self.increment_string(lastname_letters)
            name = self.base_name + str(candyear) + newname_letters
        logger.info(name)

        return name


    def is_detected_previously(self, ra, dec):
        name = xmatch_import_db(db_name=self.db_name,
                                db_user=self.db_user,
                                db_password=self.db_pwd,
                                db_table=self.cand_table_name,
                                db_output_columns=[self.db_name_field],
                                num_limit=1,
                                ra=ra,
                                dec=dec,
                                xmatch_radius_arcsec=self.xmatch_radius_arcsec,
                                db_query_columns=[],
                                db_comparison_types=[],
                                output_alias_map=[],
                                db_accepted_values=[]
                                )
        logger.info(name)
        return len(name)>0, name

    def _apply_to_candidates(
            self,
            candidate_table: pd.DataFrame,
    ) -> pd.DataFrame:
        logger.info(candidate_table)
        names = []
        lastname = None
        for ind in range(len(candidate_table)):
            cand = candidate_table.loc[ind]
            if len(cand['prv_candidates']) > 0:
                cand_name = cand['prv_candidates'][0][self.db_name_field]

            else:
                prv_det, prv_name=self.is_detected_previously(cand['ra'], cand['dec'])
                if prv_det:
                    cand_name = prv_name[0]
                else:
                    cand_name = self.get_next_name(cand[self.date_field], lastname=lastname)
                lastname = cand_name
            names.append(cand_name)
        candidate_table[self.db_name_field] = names
        logger.info(candidate_table)
        return candidate_table
