"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 24.329633712768555,
    "ZP_2.0_std": 0.13886220753192902,
    "ZP_2.0_nstars": 21,
    "ZP_3.0": 24.792879104614258,
    "ZP_3.0_std": 0.09254804253578186,
    "ZP_3.0_nstars": 21,
    "ZP_4.0": 24.95969581604004,
    "ZP_4.0_std": 0.08394473046064377,
    "ZP_4.0_nstars": 21,
    "ZP_5.0": 25.011478424072266,
    "ZP_5.0_std": 0.08213628828525543,
    "ZP_5.0_nstars": 21,
    "ZP_6.0": 25.02835464477539,
    "ZP_6.0_std": 0.08358758687973022,
    "ZP_6.0_nstars": 21,
    "ZP_7.0": 25.03775405883789,
    "ZP_7.0_std": 0.08553951978683472,
    "ZP_7.0_nstars": 21,
    "ZP_8.0": 25.0433406829834,
    "ZP_8.0_std": 0.08647112548351288,
    "ZP_8.0_nstars": 21,
    "ZP_AUTO": 25.03738784790039,
    "ZP_AUTO_std": 0.08658012002706528,
    "ZP_AUTO_nstars": 21,
    "SCORMEAN": 0.00981740837903036,
    "SCORMED": 0.07135017134681673,
    "SCORSTD": 1.2414138840965068,
}

expected_dataframe_values = {
    "magpsf": [
        16.662224170893523,
        16.94419439662437,
        17.108088175545294,
        17.17111393669992,
        15.763050183772425,
        15.834448330869588,
        15.795663417039947,
        15.76149060855563,
        15.873970432401054,
        16.063786466780208,
        15.825424225516066,
        15.72202109125955,
        15.719599862718246,
        15.8244449842761,
        15.840736419570094,
        16.034181259806044,
        15.996058656931568,
        16.254866659720413,
        16.724499296242563,
        16.726180188961727,
        16.710945493838064,
        17.536693856184392,
        17.301509351062098,
        16.383471473897316,
        16.448016329753106,
        16.103765565707818,
        16.568693202430303,
        16.611855935285234,
        16.65225803753996,
        16.258132945172555,
        17.11356111548285,
        16.912999393068723,
        17.27646497426506,
        16.918037901724972,
        16.76340612964472,
        16.38551361611458,
        16.41233808416277,
        16.772486479140603,
        16.944462281173998,
        16.331540546755427,
    ],
    "magap": [
        13.36414294414237,
        13.611801043412607,
        14.137575257128562,
        13.234456930090143,
        11.51923360327574,
        11.27831977283645,
        13.37491552738678,
        11.202726600000277,
        14.364318198377022,
        11.751389260622055,
        13.488466010104442,
        11.757313775494232,
        15.567214758708314,
        11.746227114633268,
        12.563479138722766,
        10.963807553350701,
        12.533769393754959,
        12.9398501250209,
        13.252059773574807,
        13.428161408335447,
        13.007438731672202,
        13.77082085849257,
        13.682839535644966,
        13.09641602900857,
        12.295150611994016,
        12.672502059432585,
        12.55404055334101,
        12.392014080366977,
        12.588877094142353,
        13.286369370754898,
        12.75707147152275,
        13.575407389494755,
        13.302009775222523,
        13.232048599560345,
        12.59379353443859,
        17.13041228309381,
        14.901450845973987,
        12.648822472901582,
        13.094582697943709,
        13.931117541512865,
    ],
}


pipeline = get_pipeline(
    instrument="winter", selected_configurations=["test"], night="20230726"
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

        # Expect one dataset, for one different sub-boards
        self.assertEqual(len(res[0]), 1)

        source_table = res[0][0]

        # # Uncomment to print new expected ZP dict
        print("New Results WINTER:")
        new_exp = "expected_zp = { \n"
        for header_key in source_table.keys():
            if header_key in expected_zp.keys():
                new_exp += f'    "{header_key}": {source_table[header_key]}, \n'
        new_exp += "}"
        print(new_exp)

        new_candidates_table = source_table.get_data()

        new_exp_dataframe = "expected_dataframe_values = { \n"
        for key in expected_dataframe_values:
            new_exp_dataframe += f'    "{key}": {list(new_candidates_table[key])}, \n'
        new_exp_dataframe += "}"

        print(new_exp_dataframe)

        for key, value in expected_zp.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, source_table[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, source_table[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )

        candidates_table = source_table.get_data()

        self.assertEqual(len(candidates_table), 40)
        for key, value in expected_dataframe_values.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )
