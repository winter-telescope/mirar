import logging
from glob import glob
from astropy.io import fits
from winterdrp.references.base_reference_generator import BaseReferenceGenerator
import numpy as np
from astropy.time import Time

logger = logging.getLogger(__name__)


class WIRCRef(BaseReferenceGenerator):
    abbreviation = "wirc_file_lookup"

    def __init__(
            self,
            filter_name: str,
            object_name: str,
            images_directory_path: str
    ):
        super().__init__(filter_name)
        self.filter_name = filter_name
        self.object_name = object_name
        self.images_directory_path = images_directory_path

    def get_reference(
            self,
            header: fits.Header
    ) -> fits.PrimaryHDU:
        full_imagelist = np.array(glob(f"{self.images_directory_path}/{self.object_name}/ \
                                        {self.object_name}_{self.filter_name}_*.fits"))
        imagelist = np.array([x for x in full_imagelist if not 'resamp' in x])
        try:
            mjds = [fits.getval(x, 'MJD-OBS') for x in imagelist]
        except KeyError:
            dates = Time([fits.getval(x, 'DATE') for x in imagelist])
            mjds = dates.mjd

        sinds = np.argsort(mjds)
        imagelist = imagelist[sinds]
        ref_image = imagelist[-1]

        refHDU = fits.open(ref_image)[0]
        return refHDU
