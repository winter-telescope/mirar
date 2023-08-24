"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp_list = [
    {
        "ZP_2.0": 23.89972496032715,
        "ZP_2.0_std": 0.0673225000500679,
        "ZP_2.0_nstars": 409,
        "ZP_3.0": 24.386913299560547,
        "ZP_3.0_std": 0.070705346763134,
        "ZP_3.0_nstars": 417,
        "ZP_4.0": 24.54839324951172,
        "ZP_4.0_std": 0.07586797326803207,
        "ZP_4.0_nstars": 405,
        "ZP_5.0": 24.590351104736328,
        "ZP_5.0_std": 0.08360416442155838,
        "ZP_5.0_nstars": 412,
        "ZP_6.0": 24.599721908569336,
        "ZP_6.0_std": 0.08663109689950943,
        "ZP_6.0_nstars": 405,
        "ZP_7.0": 24.603761672973633,
        "ZP_7.0_std": 0.09483309835195541,
        "ZP_7.0_nstars": 422,
        "ZP_8.0": 24.593624114990234,
        "ZP_8.0_std": 0.09740400314331055,
        "ZP_8.0_nstars": 417,
        "ZP_AUTO": 24.599775314331055,
        "ZP_AUTO_std": 0.08642666041851044,
        "ZP_AUTO_nstars": 398,
    },
    {
        "ZP_2.0": 24.176998138427734,
        "ZP_2.0_std": 0.03723874315619469,
        "ZP_2.0_nstars": 431,
        "ZP_3.0": 24.660686492919922,
        "ZP_3.0_std": 0.04019554331898689,
        "ZP_3.0_nstars": 445,
        "ZP_4.0": 24.83175277709961,
        "ZP_4.0_std": 0.040427066385746,
        "ZP_4.0_nstars": 418,
        "ZP_5.0": 24.879634857177734,
        "ZP_5.0_std": 0.04500727355480194,
        "ZP_5.0_nstars": 424,
        "ZP_6.0": 24.899553298950195,
        "ZP_6.0_std": 0.04826590046286583,
        "ZP_6.0_nstars": 427,
        "ZP_7.0": 24.911882400512695,
        "ZP_7.0_std": 0.04862583801150322,
        "ZP_7.0_nstars": 425,
        "ZP_8.0": 24.91510581970215,
        "ZP_8.0_std": 0.04960473254323006,
        "ZP_8.0_nstars": 422,
        "ZP_AUTO": 24.91059112548828,
        "ZP_AUTO_std": 0.050432268530130386,
        "ZP_AUTO_nstars": 425,
    },
]

pipeline = get_pipeline(
    instrument="winter", selected_configurations=["test"], night="20230710"
)

logging.basicConfig(level=logging.DEBUG)


class TestWinterPipeline(BaseTestCase):
    """
    Module for testing winter pipeline
    """

    def setUp(self):
        """
        Function to set up test
        Returns:

        """
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        """
        Test winter pipeline
        Returns:

        """
        self.logger.info("\n\n Testing winter pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        # Expect two datasets, for two different sub-boards
        self.assertEqual(len(res[0]), 1)
        self.assertEqual(len(res[1]), 1)

        headers = [x[0].get_header() for x in res]
        subdetids = [x["SUBDETID"] for x in headers]
        headers = [x for _, x in sorted(zip(subdetids, headers))]

        # # Uncomment to print new expected ZP dict
        print("New Results WINTER:")
        new_exp = "expected_zp_list = [ \n"
        for ind, header in enumerate(headers):
            expected_zp = expected_zp_list[ind]
            new_exp += "{ \n"
            for header_key in header.keys():
                if header_key in expected_zp.keys():
                    new_exp += f'    "{header_key}": {header[header_key]}, \n'
            new_exp += "}, \n"
        new_exp += "]"
        print(new_exp)

        for ind, header in enumerate(headers):
            expected_zp = expected_zp_list[ind]
            for key, value in expected_zp.items():
                if isinstance(value, float):
                    self.assertAlmostEqual(value, header[key], places=2)
                elif isinstance(value, int):
                    self.assertEqual(value, header[key])
                else:
                    raise TypeError(
                        f"Type for value ({type(value)} is neither float not int."
                    )
