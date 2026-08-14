from mirar.data import SourceBatch
from mirar.paths import SOURCE_NAME_KEY, TARGET_KEY
from mirar.processors.skyportal import SNCOSMO_KEY


def lmi_skyportal_formatter(source_table: SourceBatch) -> SourceBatch:
    """
    Function to add relevant fields for new sources

    :param source_table: Original source table
    :return: Updated source table
    """
    for source in source_table:
        src_df = source.get_data()

        src_df[SOURCE_NAME_KEY] = source[TARGET_KEY]
        source[SNCOSMO_KEY] = source["FILTER1"].lower().replace("-", "")

        source.set_data(src_df)

    return source_table
