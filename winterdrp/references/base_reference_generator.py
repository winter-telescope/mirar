import os
import logging
import astropy.io.fits
from winterdrp.utils.fits_tools import save_HDU_as_fits
from winterdrp.paths import base_name_key

logger = logging.getLogger(__name__)


class BaseReferenceGenerator:

    @property
    def abbreviation(self):
        raise NotImplementedError()

    def __init__(
            self,
            filter_name: str
    ):
        self.filter_name = filter_name

    @staticmethod
    def get_reference(
            header: astropy.io.fits.Header
    ) -> astropy.io.fits.PrimaryHDU:
        raise NotImplementedError()

    def write_reference(
            self,
            header: astropy.io.fits.Header,
            output_dir: str
    ) -> str:

        base_name = os.path.basename(header[base_name_key])
        logger.debug(f'Base name is {base_name}')
        refHDU = self.get_reference(
            header
        )

        output_path = self.get_output_path(output_dir, base_name).replace('.fits','') + "_ref.fits"

        # This is because Swarp requires the COADDS keyword. I am setting it to zero manually
        if not 'COADDS' in refHDU.header.keys():
            logger.debug('Setting COADDS to 0')
            refHDU.header['COADDS'] = 0
        if not 'CALSTEPS' in refHDU.header.keys():
            logger.debug('Setting CALSTEPS to blank')
            refHDU.header['CALSTEPS'] = ''

        if os.path.exists(output_path):
            os.remove(output_path)

        logger.info(f"Saving reference image to {output_path}")

        save_HDU_as_fits(refHDU, output_path)

        return output_path

    def get_output_path(
            self,
            output_dir: str,
            base_name: str
    ) -> str:
        #cat_base_name = os.path.splitext(base_name)[0] + f".{self.abbreviation}.cat"
        return os.path.join(output_dir, base_name)
