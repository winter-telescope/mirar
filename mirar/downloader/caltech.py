"""
Module for downloading data from caltech machines
"""
import logging
import os
from pathlib import Path
from typing import Optional

from mirar.paths import raw_img_dir

logger = logging.getLogger(__name__)


def download_via_ssh(
    server: str,
    base_dir: str | Path,
    night: str | int,
    pipeline: str,
    prefix: str = "",
    server_sub_dir: Optional[str] = None,
    username: str = os.getenv("SSH_USER"),
):
    """
    Function to download data via rsync+ssh from a caltech machine

    :param server: host server
    :param base_dir: base directory on host server
    :param night: night of data to download
    :param pipeline: instrument to download
    :param prefix: prefix for server data directory
    :param server_sub_dir: subdirectory within night directory on server
    :param username: username to use for ssh
    :return: None
    """
    if username is None:
        username = input(f"Please enter your username for {server}: \n")

    source_dir = f"{username}@{server}:{Path(base_dir).joinpath(prefix+night)}/"

    if server_sub_dir is not None:
        source_dir += f"{server_sub_dir}/"

    output_dir = raw_img_dir(os.path.join(pipeline, night))
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"rsync -a -v --exclude 'sub*' --exclude 'diff*' --include '*.fits' "
        f"--exclude '*' {source_dir} {output_dir}"
    )

    logger.info(f"Executing '{cmd}'")

    os.system(cmd)
