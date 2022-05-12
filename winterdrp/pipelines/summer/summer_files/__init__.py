import os
from winterdrp.pipelines.summer.summer_files.schema import get_summer_schema_path
summer_dir = os.path.dirname(__file__)

astromatic_config_dir = os.path.join(summer_dir, 'config')
summer_mask_path = os.path.join(summer_dir, "mask.fits")
summer_weight_path = os.path.join(summer_dir, "weight.fits")

sextractor_astrometry_config = {
    "config_path": os.path.join(astromatic_config_dir, 'astrom.sex'),
    "filter_path": os.path.join(astromatic_config_dir, 'default.conv'),
    "parameter_path": os.path.join(astromatic_config_dir, 'astrom.param'),
    "starnnw_path": os.path.join(astromatic_config_dir, 'default.nnw')
}


sextractor_photometry_config = {
    "config_path": os.path.join(astromatic_config_dir, 'photomCat.sex'),
    "filter_path": os.path.join(astromatic_config_dir, 'default.conv'),
    "parameter_path": os.path.join(astromatic_config_dir, 'photom.param'),
    "starnnw_path": os.path.join(astromatic_config_dir, 'default.nnw')
}

astrom_scamp = 'config/scamp.conf'
astrom_sex = 'config/astrom.sex'
astrom_param = 'config/astrom.param'
astrom_filter = 'config/default.conv'
astrom_swarp = 'config/config.swarp'
astrom_nnw = 'config/default.nnw'

scamp_path = os.path.join(astromatic_config_dir, "scamp.conf")

swarp_path = os.path.join(astromatic_config_dir, "config.swarp")
