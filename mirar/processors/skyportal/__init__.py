"""
Central module for the skyportal integration
"""

from mirar.processors.skyportal.client import SkyportalClient
from mirar.processors.skyportal.skyportal_candidate import SkyportalCandidateUploader
from mirar.processors.skyportal.skyportal_source import (
    SNCOSMO_KEY,
    SkyportalSourceUploader,
)
from mirar.processors.skyportal.thumbnail import make_thumbnail
