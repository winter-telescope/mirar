from mirar.paths import BASE_NAME_KEY
from mirar.pipelines.cfht.config import sextractor_photometry_config
from mirar.pipelines.cfht.generator import cfht_anet_sextractor_config_path_generator
from mirar.pipelines.cfht.load_cfht_file import load_raw_cfht_image
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.astrometry.anet import AstrometryNet
from mirar.processors.utils import (
    CustomImageBatchModifier,
    HeaderAnnotator,
    HeaderEditor,
    ImageBatcher,
    ImageDebatcher,
    ImageLoader,
    ImagePlotter,
    ImageRebatcher,
    ImageRejector,
    ImageSaver,
    ImageSelector,
    MEFLoader,
)

astrometry = [
    ImageLoader(input_sub_dir="single_ext", load_image=load_raw_cfht_image),
    ImageBatcher(BASE_NAME_KEY),
    AstrometryNet(
        output_sub_dir="anet",
        scale_bounds=[0.1, 0.3],
        scale_units="app",
        use_sextractor=False,
        timeout=120,
        search_radius_deg=1.0,
        cache=False,
        x_image_key=None,
        y_image_key=None,
        downsample=2,
        write_regions=False,
    ),
    ImageSaver(output_dir_name="post_anet"),
    Sextractor(
        **sextractor_photometry_config,
        output_sub_dir="photometery",
        write_regions_bool=True,
    ),
]
