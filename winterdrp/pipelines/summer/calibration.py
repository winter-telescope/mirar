import pandas as pd


def select_bias(
       observing_log: pd.DataFrame
) -> [str]:
    mask = observing_log["OBSTYPE"] == "BIAS"
    return list(observing_log[mask]["RAWIMAGEPATH"])


def select_flats_archival(
       observing_log: pd.DataFrame
) -> [str]:
    mask = observing_log["OBSTYPE"] == "FLAT"
    return list(observing_log[mask]["RAWIMAGEPATH"])
