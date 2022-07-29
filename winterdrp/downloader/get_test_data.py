import wget
import os
from winterdrp.paths import raw_img_dir

TEST_DATA_URL = "https://raw.githubusercontent.com/robertdstein/wircpipe/raw/master/master_flats/master_flat_H.fits"


def download_test_data(
        url: str = TEST_DATA_URL
):

    test_data_dir = raw_img_dir("summer/testdata/")

    try:
        os.makedirs(test_data_dir)
    except OSError:
        pass

    output_path = os.path.join(test_data_dir, os.path.basename(TEST_DATA_URL))

    wget.download(TEST_DATA_URL, output_path)


download_test_data()
