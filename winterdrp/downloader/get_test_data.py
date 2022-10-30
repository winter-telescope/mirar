import os
from winterdrp.paths import winter_code_dir
import logging

logger = logging.getLogger(__name__)

TEST_DATA_URL = "git@github.com:winter-telescope/wirc_starterpack.git"

test_data_dir = os.path.join(
    os.path.dirname(winter_code_dir),
    os.path.basename(TEST_DATA_URL.replace(".git", ""))
)

TEST_DATA_TAG = "v0.1.3"


def get_test_data_dir() -> str:

    if not os.path.isdir(test_data_dir):

        cmd = f"git clone {TEST_DATA_URL} {test_data_dir}"

        logger.info(f"No test data found. Downloading. Executing: {cmd}")

        os.system(cmd)

    else:
        cmds = [
            f"git -C {test_data_dir} checkout main",
            f"git -C {test_data_dir} pull"
        ]

        for cmd in cmds:
            logger.info(f"Trying to update test data. Executing: {cmd}")
            os.system(cmd)

    fix_version_cmd = f"git -C {test_data_dir} checkout -d tags/{TEST_DATA_TAG}"

    logger.info(f"Checkout out correct test data commit. Executing: {fix_version_cmd}")

    os.system(fix_version_cmd)

    return test_data_dir


if __name__ == "__main__":
    get_test_data_dir()
