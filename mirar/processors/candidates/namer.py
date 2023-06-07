"""
Module containing a processor for assigning names to sources
"""
import logging

from astropy.time import Time

from mirar.data import SourceBatch
from mirar.paths import CAND_NAME_KEY
from mirar.processors.base_processor import BaseDataframeProcessor
from mirar.processors.database import BaseDatabaseProcessor

logger = logging.getLogger(__name__)


class CandidateNamer(BaseDatabaseProcessor, BaseDataframeProcessor):
    """Processor to sequentially assign names to sources, of the form a, aa, aba..."""

    base_key = "namer"

    # Go one at a time to avoid... race conditions
    max_n_cpu = 1

    def __init__(
        self,
        base_name: str,
        xmatch_radius_arcsec: float,
        name_start: str = "aaaaa",
        db_name_field: str = CAND_NAME_KEY,
        db_order_field: str = "candid",
        date_field: str = "jd",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.db_name_field = db_name_field
        self.db_order_field = db_order_field
        self.base_name = base_name
        self.name_start = name_start
        self.date_field = date_field
        self.crossmatch_radius_arcsec = xmatch_radius_arcsec

    @staticmethod
    def increment_string(string: str):
        """

        Parameters
        ----------
        string

        Returns
        -------
        An incremented string, eg. aaa -> aab, aaz -> aba, azz -> baa, zzz-> aaaa
        """
        character_position = len(string) - 1
        # will iteratively try to increment characters starting from the last
        increment_bool = False
        new_string = ""
        while character_position >= 0:
            cref = string[character_position]
            if increment_bool:
                new_string = cref + new_string
                character_position -= 1
                continue
            cref_ordered = ord(cref)
            # increment each character, if at 'z', increment the next one
            if cref_ordered + 1 > 122:
                new_string = "a" + new_string
                if character_position == 0:
                    new_string = "a" + new_string
            else:
                next_character = chr(cref_ordered + 1)
                new_string = next_character + new_string
                increment_bool = True
            character_position -= 1
            continue

        return new_string

    def get_next_name(self, cand_jd: float, last_name: str = None) -> str:
        """
        Function to get a new candidate name

        :param cand_jd: jd of detection
        :param last_name: last name
        :return: new name
        """
        cand_year = Time(cand_jd, format="jd").datetime.year % 1000
        if last_name is None:
            query = (
                f'SELECT "{self.db_name_field}" FROM {self.db_table} '
                f"ORDER BY {self.db_order_field} desc LIMIT 1;"
            )

            res = self.pg_user.execute_query(query, db_name=self.db_name)

            if len(res) == 0:
                name = self.base_name + str(cand_year) + self.name_start
                return name

            last_name = res[0][0]
            logger.debug(res)

        last_year = int(last_name[len(self.base_name) : len(self.base_name) + 2])
        if cand_year != last_year:
            name = self.base_name + str(cand_year) + self.name_start
        else:
            last_name_letters = last_name[len(self.base_name) + 2 :]
            new_name_letters = self.increment_string(last_name_letters)
            name = self.base_name + str(cand_year) + new_name_letters
        logger.debug(name)
        return name

    def is_detected_previously(self, ra_deg: float, dec_deg: float) -> tuple[bool, str]:
        """
        Checks whether a source has been detected previously

        :param ra_deg: ra (deg)
        :param dec_deg: dec (deg)
        :return: boolean whether a source has been detected previously
        """
        name = self.pg_user.crossmatch_with_database(
            db_name=self.db_name,
            db_table=self.db_table,
            db_output_columns=[self.db_name_field],
            num_limit=1,
            ra=ra_deg,
            dec=dec_deg,
            crossmatch_radius_arcsec=self.crossmatch_radius_arcsec,
        )  # [0]
        logger.info(name)
        return len(name) > 0, name

    def _apply_to_candidates(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        for source_table in batch:
            candidate_table = source_table.get_data()

            names = []
            lastname = None
            for ind in range(len(candidate_table)):
                cand = candidate_table.loc[ind]
                if len(cand["prv_candidates"]) > 0:
                    cand_name = cand["prv_candidates"][0][self.db_name_field]

                else:
                    prv_det, prv_name = self.is_detected_previously(
                        cand["ra"], cand["dec"]
                    )

                    if prv_det:
                        cand_name = prv_name
                    else:
                        cand_name = self.get_next_name(
                            cand[self.date_field], last_name=lastname
                        )
                    lastname = cand_name
                names.append(cand_name)

            candidate_table[self.db_name_field] = names
            source_table.set_data(candidate_table)

        return batch
