import logging
import os
import logging
from winterdrp.processors.base_processor import BaseDataframeProcessor
import numpy as np
import pandas as pd
from winterdrp.processors.database.postgres import execute_query, xmatch_import_db
from astropy.time import Time
from time import sleep

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
                 naming_sequence_name: str = 'name_count',
                 year_sequence_name: str = 'year_count',
                 num_letters : int = 5,
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
        self.naming_sequence_name = naming_sequence_name
        self.year_sequence_name = year_sequence_name
        self.num_letters = num_letters

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

    @staticmethod
    def numberToBase(n, b):
        if n == 0:
            return [0]
        digits = []
        while n:
            digits.append(int(n % b))
            n //= b
        return digits[::-1]

    def get_char_name_from_num(self, number):
        base27 = self.numberToBase(number-1, 27)
        text = [chr(97 + x) for x in base27]
        name = "".join(text).rjust(self.num_letters, 'a')
        return name

    def get_next_num_year(self):
        query = f"""SELECT nextval('{self.naming_sequence_name}');"""
        res = execute_query(query, db_name=self.db_name, db_user=self.db_user, password=self.db_pwd)
        nextseq_num = res[0][0]

        query = f"""SELECT last_value FROM {self.year_sequence_name};"""
        res = execute_query(query, db_name=self.db_name, db_user=self.db_user, password=self.db_pwd)
        naming_cur_year = res[0][0]
        return nextseq_num, naming_cur_year

    def get_next_name(self, candjd, lastname=None):
        candyear = Time(candjd, format='jd').datetime.year % 100

        nextseq_num, naming_cur_year = self.get_next_num_year()

        assert naming_cur_year <= candyear

        if candyear != naming_cur_year:
            niter = 5
            year_incremented = False
            for i in range(niter):
                # This needs to be supplemented with a cron that runs every day at noon PT to update the year counter,
                # to avoid any potential race conditions

                # Wait for 25 seconds in case some other process resets the year
                requery_nextseq_num, requery_naming_cur_year = self.get_next_num_year()
                if requery_naming_cur_year == candyear:
                    nextseq_num = requery_nextseq_num
                    naming_cur_year = requery_naming_cur_year
                    year_incremented = True
                    break
                elif requery_naming_cur_year > candyear:
                    query = f"""ALTER SEQUENCE {self.year_sequence_name} RESTART {candyear};"""
                    res = execute_query(query, db_name=self.db_name, db_user=self.db_user, password=self.db_pwd)
                    year_incremented = True
                    break
                sleep(5)

            if not year_incremented:
                query = f"""ALTER SEQUENCE {self.naming_sequence_name} RESTART {1};"""
                res = execute_query(query, db_name=self.db_name, db_user=self.db_user, password=self.db_pwd)

                query = f"""SELECT nextval('{self.year_sequence_name}');"""
                res = execute_query(query, db_name=self.db_name, db_user=self.db_user, password=self.db_pwd)

                nextseq_num, naming_cur_year = self.get_next_num_year()
                assert naming_cur_year == candyear

        next_char_name = self.get_char_name_from_num(int(str(nextseq_num)[2:]))
        name = str(naming_cur_year)+next_char_name

        '''
        if lastname is None:
            # query = f"""SELECT "{self.db_name_field}" FROM {self.cand_table_name} ORDER BY {self.db_order_field} desc LIMIT 1;"""
            query = f"""SELECT nextval('{self.naming_sequence_name}');"""
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
        '''
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
        return len(name) > 0, name

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
                prv_det, prv_name = self.is_detected_previously(cand['ra'], cand['dec'])
                if prv_det:
                    cand_name = prv_name[0]
                else:
                    cand_name = self.get_next_name(cand[self.date_field], lastname=lastname)
                lastname = cand_name
            names.append(cand_name)
        candidate_table[self.db_name_field] = names
        logger.info(candidate_table)
        return candidate_table
