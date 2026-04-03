"""
Constants for the winter pipeline.
"""

import astroplan
import astropy.coordinates as coords

mirage_filters_map = {"Y": 1, "J": 2, "Hs": 3, "dark": 4}

mirage_inv_filters_map = {v: k for k, v in mirage_filters_map.items()}

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

PALOMAR_LOC = coords.EarthLocation(
    lat=coords.Latitude("33d21m25.5s"),
    lon=coords.Longitude("-116d51m58.4s"),
    height=1696.0,
)
palomar_observer = astroplan.Observer(location=PALOMAR_LOC)
