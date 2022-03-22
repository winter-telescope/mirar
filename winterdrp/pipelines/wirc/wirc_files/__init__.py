import os
from winterdrp.processors.astromatic.config import astromatic_config_dir

wirc_dir = os.path.dirname(__file__)

wirc_mask_path = os.path.join(wirc_dir, "wirc_bad_feb_2018.fits")

sextractor_astrometry_config = {
    "config_path": os.path.join(wirc_dir, 'stack.sex'),
    "filter_path": os.path.join(astromatic_config_dir, 'default.conv'),
    "parameter_path": os.path.join(astromatic_config_dir, 'default.param'),
    "starnnw_path": os.path.join(astromatic_config_dir, 'default.nnw')
}