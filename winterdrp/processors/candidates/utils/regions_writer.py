import os
import pandas as pd
from winterdrp.processors.base_processor import BaseDataframeProcessor
from winterdrp.paths import base_output_dir, get_output_dir, get_output_path
import logging

logger = logging.getLogger(__name__)


class RegionsWriter(BaseDataframeProcessor):

    def __init__(self,
                 output_dir_name: str = None,
                 region_pix_radius: float = 8,
                 output_dir: str = base_output_dir,
                 *args,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.output_dir_name = output_dir_name
        self.region_pix_radius = region_pix_radius
        self.output_dir = output_dir

    def _apply_to_candidates(
            self,
            candidate_table: pd.DataFrame,
    ) -> pd.DataFrame:
        try:
            os.makedirs(get_output_dir(
                dir_root=self.output_dir_name,
                sub_dir=self.night_sub_dir,
                output_dir=self.output_dir
            ))
        except OSError:
            pass
        regions_basepath = os.path.basename(candidate_table.iloc[0]['diffimname']).replace('.fits', '.reg')
        regions_path = get_output_path(regions_basepath,
                                       dir_root=self.output_dir_name,
                                       sub_dir=self.night_sub_dir,
                                       output_dir=self.output_dir
                                       )
        logger.info(f'Writing regions path to {regions_path}')
        with open(f"{regions_path}", 'w') as f:
            f.write('image\n')
            for ind in range(len(candidate_table)):
                row = candidate_table.iloc[ind]
                f.write(f"CIRCLE({row['X_IMAGE']},{row['Y_IMAGE']},{self.region_pix_radius})\n")

        return candidate_table
