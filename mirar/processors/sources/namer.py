"""
Module containing a processor for assigning names to sources
"""

import logging

import pandas as pd
from astropy.time import Time
from sqlalchemy import select

from mirar.data import SourceBatch
from mirar.database.transactions.select import run_select
from mirar.paths import SOURCE_NAME_KEY, TIME_KEY
from mirar.processors.database.database_selector import SingleSpatialCrossmatchSource

logger = logging.getLogger(__name__)


class CandidateNamer(SingleSpatialCrossmatchSource):
    """Processor to sequentially assign names to sources, of the form a, aa, aba..."""

    base_key = "namer"

    # Go one at a time to avoid... race conditions
    max_n_cpu = 1

    def __init__(
        self,
        base_name: str,
        name_start: str = "aaaaa",
        db_name_field: str = SOURCE_NAME_KEY,
        name_key: str = SOURCE_NAME_KEY,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.db_name_field = db_name_field

        # Ensure that the name field is in the output columns
        self.db_output_columns = list(
            set([self.db_name_field] + self.db_output_columns)
        )

        self.base_name = base_name
        self.name_start = name_start
        self.name_key = name_key
        self.lastname = None

    def description(self) -> str:
        return (
            f"Sequentially assign names to new sources, e.g "
            f"{self.base_name}24{self.name_start}"
        )

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

    def extract_last_year(self, last_name: str) -> int:
        """
        Extract the year from the last name

        :param last_name: last name
        :return: year
        """
        last_year = int(last_name[len(self.base_name) : len(self.base_name) + 2])
        return last_year

    def get_next_name(self, detection_time: Time, last_name: str = None) -> str:
        """
        Function to get a new candidate name

        :param detection_time: detection time (Astropy Time object)
        :param last_name: last name
        :return: new name
        """
        cand_year = detection_time.datetime.year % 1000

        if last_name is not None:
            last_year = self.extract_last_year(last_name)
            if last_year != cand_year:
                last_name = None

        if last_name is None:

            col = self.db_table.sql_model.__table__.c[self.db_name_field]

            # Select most recent name of same year
            sel = (
                select(col).where(col.contains(cand_year)).order_by(col.desc()).limit(1)
            )

            res = run_select(query=sel, sql_table=self.db_table.sql_model)

            # If no names of the same year, start from the beginning
            if len(res) == 0:
                name = self.base_name + str(cand_year) + self.name_start
                return name

            last_name = res[self.db_name_field].iloc[0]
            logger.debug(res)

        last_year = self.extract_last_year(last_name)

        assert (
            last_year == cand_year
        ), f"Last year {last_year} does not match candidate year {cand_year}"

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

            metadata = source_table.get_metadata()

            detection_time = Time(source_table[TIME_KEY])

            matches = []

            for ind, source in sources.iterrows():

                match = self.query_for_source(source, metadata)

                if len(match) > 0:
                    source_name = match[self.name_key].iloc[0]
                    logger.debug(f"Source already has name: {source_name}")
                    matches.append(match)

                else:
                    source_name = self.get_next_name(
                        detection_time, last_name=self.lastname
                    )
                    self.lastname = source_name
                    logger.debug(f"Assigning name: {source_name} to source # {ind}.")

                    # Insert the name into the source table
                    source[self.name_key] = source_name

                    metadata_dict = self.generate_super_dict(metadata, source)

                    new = self.db_table(**metadata_dict)
                    match = new.insert_entry(
                        duplicate_protocol="fail",
                        returning_key_names=self.db_output_columns,
                    )
                    matches.append(match)

            match_df = pd.concat(matches, ignore_index=True, axis=0)

            for column in self.db_output_columns:
                sources[column] = match_df[column]

            source_table.set_data(sources)

        return batch
