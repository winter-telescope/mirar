"""
Constants for the winter pipeline.
"""
import astroplan
import astropy.coordinates as coords
import pandas as pd

winter_filters_map = {"Y": 1, "J": 2, "Hs": 3, "dark": 4}

imgtype_dict = {
    "science": "SCIENCE",
    "bias": "CAL",
    "flat": "CAL",
    "dark": "CAL",
    "focus": "FOCUS",
    "pointing": "POINTING",
    "other": "NULL",
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

PALOMAR_LOC = coords.EarthLocation(
    lat=coords.Latitude("33d21m25.5s"),
    lon=coords.Longitude("-116d51m58.4s"),
    height=1696.0,
)
palomar_observer = astroplan.Observer(location=PALOMAR_LOC)
