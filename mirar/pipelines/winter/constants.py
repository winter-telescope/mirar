"""
Constants for the winter pipeline.
"""
import pandas as pd

winter_filters_map = {"Y": 1, "J": 2, "Hs": 3, "dark": 4}

imgtype_dict = {
    "SCIENCE": "SCIENCE",
    "BIAS": "CAL",
    "FLAT": "CAL",
    "DARK": "CAL",
    "FOCUS": "FOCUS",
    "POINTING": "POINTING",
    "OTHER": "NULL",
}

NXSPLIT = 1
NYSPLIT = 2

WINTER_N_BOARDS = 6

_subdets = []

for ndetector in range(WINTER_N_BOARDS):
    for nx in range(NXSPLIT):
        for ny in range(NYSPLIT):
            _subdets.append(
                {
                    "boardid": ndetector,
                    "n_board_max": WINTER_N_BOARDS,
                    "nx": nx + 1,
                    "nxtot": NXSPLIT,
                    "ny": ny + 1,
                    "nytot": NYSPLIT,
                }
            )

subdets = pd.DataFrame(_subdets)
subdets["subdetid"] = range(1, len(subdets) + 1)
