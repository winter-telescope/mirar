import os
from winterdrp.paths import astrometry_output_dir
from winterdrp.calibrate.sourceextractor import run_sextractor
from winterdrp.preprocessing.base_processor import BaseProcessor


class SextractorRunner(BaseProcessor):

    base_key = "sextractor"

    requires = ["save"]

    def _apply_to_images(
            self,
            images: list,
            headers: list,
            sub_dir: str = ""
    ) -> (list, list):

        # Try making output directory, unless it exists

        output_dir = astrometry_output_dir(sub_dir)

        try:
            os.makedirs(output_dir)
        except OSError:
            pass

        for i, image in enumerate(images):

            # First run Sextractor

            run_sextractor(
                image,
                output_dir=output_dir,
            )

        return images, headers
