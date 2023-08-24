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
        "ZP_2.0": 23.90024757385254,
        "ZP_2.0_std": 0.06780964881181717,
        "ZP_2.0_nstars": 411,
        "ZP_3.0": 24.386281967163086,
        "ZP_3.0_std": 0.0703275129199028,
        "ZP_3.0_nstars": 415,
        "ZP_4.0": 24.547704696655273,
        "ZP_4.0_std": 0.07585957646369934,
        "ZP_4.0_nstars": 405,
        "ZP_5.0": 24.589977264404297,
        "ZP_5.0_std": 0.0833786204457283,
        "ZP_5.0_nstars": 411,
        "ZP_6.0": 24.60162353515625,
        "ZP_6.0_std": 0.08639319241046906,
        "ZP_6.0_nstars": 404,
        "ZP_7.0": 24.604106903076172,
        "ZP_7.0_std": 0.09461381286382675,
        "ZP_7.0_nstars": 421,
        "ZP_8.0": 24.594432830810547,
        "ZP_8.0_std": 0.09780757129192352,
        "ZP_8.0_nstars": 419,
        "ZP_AUTO": 24.599721908569336,
        "ZP_AUTO_std": 0.08646132051944733,
        "ZP_AUTO_nstars": 398,
    },
    {
        "ZP_2.0": 24.1755313873291,
        "ZP_2.0_std": 0.03563622757792473,
        "ZP_2.0_nstars": 422,
        "ZP_3.0": 24.661762237548828,
        "ZP_3.0_std": 0.03949178382754326,
        "ZP_3.0_nstars": 436,
        "ZP_4.0": 24.835601806640625,
        "ZP_4.0_std": 0.04035663977265358,
        "ZP_4.0_nstars": 414,
        "ZP_5.0": 24.880395889282227,
        "ZP_5.0_std": 0.04603949934244156,
        "ZP_5.0_nstars": 427,
        "ZP_6.0": 24.89640235900879,
        "ZP_6.0_std": 0.04970197379589081,
        "ZP_6.0_nstars": 431,
        "ZP_7.0": 24.910377502441406,
        "ZP_7.0_std": 0.05080369487404823,
        "ZP_7.0_nstars": 432,
        "ZP_8.0": 24.91377067565918,
        "ZP_8.0_std": 0.04984834045171738,
        "ZP_8.0_nstars": 419,
        "ZP_AUTO": 24.909290313720703,
        "ZP_AUTO_std": 0.05263081565499306,
        "ZP_AUTO_nstars": 433,
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
