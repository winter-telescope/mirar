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

def check_gitlfs():
    try:
        execute("which git-lfs")
    except ExecutionError:
        raise FileNotFoundError("Downloading test data requires git-lfs. "
                                "See install instructions: https://git-lfs.github.com/")


def get_wirc_test_data() -> str:

    if not os.path.isdir(test_data_dir):

        cmd = f"git clone {TEST_DATA_URL} {test_data_dir}"

        logger.info(f"No test data found. Downloading. Executing: {cmd}")

        check_gitlfs()

        os.system(cmd)

    zip_data_path = glob(f"{test_data_dir}/*.zip")
    for path in zip_data_path:
        logger.info(f"Unzipping {path}")
        with zipfile.ZipFile(path, 'r') as zip_ref:
            zip_ref.extractall(test_data_dir)

    sub_dirs = [
        x for x in glob(f"{test_data_dir}/*")
        if np.logical_and(os.path.isdir(x), os.path.basename(x)[0] not in [".", "_"])
    ]
    assert len(sub_dirs) == 1
    test_data_path = sub_dirs[0]

    return test_data_path
