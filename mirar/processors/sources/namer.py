"""
Module containing a processor for assigning names to sources
"""

import logging

import numpy as np
from astropy.time import Time
from sqlalchemy import select, text

from mirar.data import SourceBatch
from mirar.database.transactions.select import run_select
from mirar.paths import SOURCE_NAME_KEY, SOURCE_XMATCH_KEY, TIME_KEY
from mirar.processors.base_processor import PrerequisiteError
from mirar.processors.database import CrossmatchSourceWithDatabase
from mirar.processors.database.database_selector import BaseDatabaseSourceSelector

logger = logging.getLogger(__name__)


class CandidateNamer(BaseDatabaseSourceSelector):
    """Processor to sequentially assign names to sources, of the form a, aa, aba..."""

    base_key = "namer"

    # Go one at a time to avoid... race conditions
    max_n_cpu = 1

    def __init__(
        self,
        base_name: str,
        name_start: str = "aaaaa",
        db_name_field: str = SOURCE_NAME_KEY,
        db_order_field: str = "candid",
        **kwargs,
    ):
        super().__init__(db_output_columns=[db_name_field], **kwargs)
        self.db_name_field = db_name_field
        self.db_order_field = db_order_field
        self.base_name = base_name
        self.name_start = name_start
        self.lastname = None

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

    def get_next_name(self, detection_time: Time, last_name: str = None) -> str:
        """
        Function to get a new candidate name

        :param detection_time: detection time (Astropy Time object)
        :param last_name: last name
        :return: new name
        """
        cand_year = detection_time.datetime.year % 1000
        if last_name is None:
            res = run_select(
                query=select(getattr(self.db_table.sql_model, self.db_name_field))
                .order_by(text(f"{self.db_order_field} desc"))
                .limit(1),
                sql_table=self.db_table.sql_model,
            )

            if len(res) == 0:
                name = self.base_name + str(cand_year) + self.name_start
                return name

            last_name = res[SOURCE_NAME_KEY].iloc[0]
            logger.debug(res)

        last_year = int(last_name[len(self.base_name) : len(self.base_name) + 2])
        if cand_year != last_year:
            name = self.base_name + str(cand_year) + self.name_start
        else:
            last_name_letters = last_name[len(self.base_name) + 2 :]
            new_name_letters = self.increment_string(last_name_letters)
            name = self.base_name + str(cand_year) + new_name_letters
        logger.debug(f"Assigning name: {name}")
        return name

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        for source_table in batch:
            sources = source_table.get_data()

            assert (
                SOURCE_XMATCH_KEY in sources.columns
            ), "No candidate cross-match in source table"

            names = []

            detection_time = Time(source_table[TIME_KEY])
            for ind, source in sources.iterrows():
                if len(source[SOURCE_XMATCH_KEY]) > 0:
                    source_name = source[SOURCE_XMATCH_KEY][0][self.db_name_field]
                else:
                    source_name = self.get_next_name(
                        detection_time, last_name=self.lastname
                    )
                    self.lastname = source_name
                logger.debug(f"Assigning name: {source_name} to source # {ind}.")
                names.append(source_name)

            sources[self.db_name_field] = names
            source_table.set_data(sources)

        return batch

    def check_prerequisites(
        self,
    ):
        check = np.sum(
            [isinstance(x, CrossmatchSourceWithDatabase) for x in self.preceding_steps]
        )
        if check < 1:
            err = (
                f"{self.__module__} requires {CrossmatchSourceWithDatabase} "
                f"as a prerequisite. "
                f"However, the following steps were found: {self.preceding_steps}."
            )
            logger.error(err)
            raise PrerequisiteError(err)
