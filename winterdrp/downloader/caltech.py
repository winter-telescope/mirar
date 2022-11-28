import os
import logging
import os
from pathlib import Path
from winterdrp.paths import raw_img_dir


logger = logging.getLogger(__name__)


def download_via_ssh(
        server: str,
        base_dir: str | Path,
        night: str | int,
        pipeline: str,
        prefix: str = "",
        server_sub_dir: str = None,
        username: str = os.getenv("SSH_USER")
):
    if username is None:
        username = input(f"Please enter your username for {server}: \n")

    source_dir = f"{username}@{server}:{Path(base_dir).joinpath(prefix+night)}/"

    if server_sub_dir is not None:
        source_dir += f"{server_sub_dir}/"

    output_dir = raw_img_dir(os.path.join(pipeline, night))

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    cmd = f"rsync -a -v --exclude 'sub*' --exclude 'diff*' --include '*.fits' --exclude '*' {source_dir} {output_dir}"

    logger.info(f"Executing '{cmd}'")

    os.system(cmd)

