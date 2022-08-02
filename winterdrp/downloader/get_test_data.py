import wget
import os
from winterdrp.paths import winter_code_dir
import logging
from glob import glob
import zipfile
import numpy as np
from winterdrp.utils import execute, ExecutionError

logger = logging.getLogger(__name__)

TEST_DATA_URL = "git@github.com:winter-telescope/wirc_starterpack.git"

test_data_dir = os.path.join(
    os.path.dirname(winter_code_dir),
    os.path.basename(TEST_DATA_URL.replace(".git", ""))
)


def get_test_data() -> str:

    if not os.path.isdir(test_data_dir):

        cmd = f"git clone {TEST_DATA_URL} {test_data_dir}"

        logger.info(f"No test data found. Downloading. Executing: {cmd}")

        os.system(cmd)

    return test_data_dir


if __name__ == "__main__":
    get_test_data()
