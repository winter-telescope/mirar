"""
Constants for the winter pipeline.
"""

import astroplan
import astropy.coordinates as coords
import pandas as pd

winter_filters_map = {"Y": 1, "J": 2, "Hs": 3, "dark": 4}

winter_inv_filters_map = {v: k for k, v in winter_filters_map.items()}

sncosmo_filters = {
    "y": "desy",
    "j": "2massj",
    "h": "2massh",
}

imgtype_dict = {
    "science": "SCIENCE",
    "bias": "CAL",
    "flat": "CAL",
    "dark": "CAL",
    "focus": "FOCUS",
    "pointing": "POINTING",
    "other": "NULL",
    "corrupted": "CORRUPTED",
}

# Typical radial offset of a single WINTER dither from the nominal field
# center. Used as a safety margin when sizing reference-catalog queries, so
# a single query/cache (keyed by field/subdet/filter, shared across all
# dithers of one visit) has enough coverage for every dither, not just the
# one it happened to be centered on. Empirically, offsets in real WINTER
# dither patterns run up to ~3 arcmin from center (see e.g. the 20230726
# test dataset); this adds a healthy margin above that.
WINTER_DITHER_RADIUS_ARCMIN = 5.0

NXSPLIT = 1
NYSPLIT = 1

WINTER_N_BOARDS = 6

_subdets = []

all_winter_board_ids = [0, 1, 2, 3, 4, 5, 6]

for ndetector in all_winter_board_ids:
    for nx in range(NXSPLIT):
        for ny in range(NYSPLIT):
            _subdets.append(
                {
                    "boardid": ndetector,
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
