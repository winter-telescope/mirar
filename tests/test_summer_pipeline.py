import unittest
import logging
from winterdrp.pipelines import get_pipeline
from winterdrp.data import Dataset, ImageBatch

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 24.381802801640834,
    "ZP_2.0_std": 0.08222868531222044,
    "ZP_2.0_nstars": 30,
    "ZP_3.0": 25.06422489916484,
    "ZP_3.0_std": 0.06870406816807186,
    "ZP_3.0_nstars": 30,
    "ZP_4.0": 25.451230355453497,
    "ZP_4.0_std": 0.06158930864995796,
    "ZP_4.0_nstars": 30,
    "ZP_5.0": 25.659666622034713,
    "ZP_5.0_std": 0.06283258570749543,
    "ZP_5.0_nstars": 30,
    "ZP_6.0": 25.78091181551616,
    "ZP_6.0_std": 0.06342292306060344,
    "ZP_6.0_nstars": 30,
    "ZP_7.0": 25.853334717305508,
    "ZP_7.0_std": 0.06322232852162107,
    "ZP_7.0_nstars": 30,
    "ZP_8.0": 25.90134099756877,
    "ZP_8.0_std": 0.06359669923996372,
    "ZP_8.0_nstars": 30
}

test_config_name = "test"

pipeline = get_pipeline(
    instrument="summer",
    selected_configurations=["test"],
    night="20220402"
)


class TestSummerPipeline(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        self.logger.info("\n\n Testing summer pipeline \n\n")

        res, errorstack = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        self.assertEqual(len(res[0]), 1)

        header = res[0][0].get_header()

        for key, value in expected_zp.items():
            if isinstance(value, float):
                print(key, value, header[key])
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(f"Type for value ({type(value)} is neither float not int.")


