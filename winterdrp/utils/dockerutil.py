"""
Module containing docker integration (beta-stage)
"""
import io
import logging
import os
import tarfile
from pathlib import Path

import docker
from docker.errors import DockerException
from docker.models.containers import Container

logger = logging.getLogger(__name__)

DOCKER_IMAGE_NAME = "robertdstein/astrodocker"
docker_dir = Path("/usr/src/astrodocker")


def new_container():
    """Generate a new docker.models.containers.Container object, using the default
    docker daemon and the docker image.
    If the image is not found locally,
    the image will first be pulled from DockerHub.

    This function requires a Docker daemon to first be running.

    Returns
    -------
    A docker container built with the {DOCKER_IMAGE_NAME} image
    """
    try:
        client = docker.from_env()
    except DockerException as exc:
        err = (
            "Unable to connect to Docker daemon. "
            "Have you installed Docker, and started a daemon? "
            "Find out more at https://www.docker.com "
        )
        logger.error(err)
        raise ConnectionError(err) from exc

    if len(client.images.list(DOCKER_IMAGE_NAME)) < 1:
        logger.info(f"Pulling docker image {DOCKER_IMAGE_NAME}")
        client.images.pull(DOCKER_IMAGE_NAME)

    return client.containers.run(DOCKER_IMAGE_NAME, tty=True, detach=True)


def docker_path(file_path: str | Path) -> Path:
    """
    Converts a local path to the corresponding path in the docker container

    :param file_path: file path
    :return:
    """
    return docker_dir.joinpath(Path(file_path).name)


def docker_get(container: Container, local_path: str | Path):
    """Function to cope one file from the Docker container 'container' to 'local_path'.
    The file in the container should have
    the same name as the base file in 'local_path'.

    Parameters
    ----------
    container: A docker.models.container.Container object
    local_path: Local path of file to copy to

    Returns
    -------
    """

    container_path = docker_path(local_path)

    with open(local_path, "wb") as local_file:
        bits, _ = container.get_archive(container_path.as_posix())
        for chunk in bits:
            local_file.write(chunk)


def docker_put(container: Container, local_path: str | Path):
    """Function to one file, at 'local_path' into the Docker container 'container'

    Parameters
    ----------
    container: A docker.models.container.Container object
    local_path: Local path of file to copy

    Returns
    -------
    """
    stream = io.BytesIO()

    with tarfile.open(fileobj=stream, mode="w|") as tar, open(
        local_path, "rb"
    ) as local_file:
        info = tar.gettarinfo(fileobj=local_file)
        info.name = os.path.basename(local_path)
        tar.addfile(info, local_file)

    return container.put_archive(docker_dir.as_posix(), stream.getvalue())


def docker_batch_put(container: Container, local_paths: str | list):
    """Function to copy multiple files into a Docker container

    Parameters
    ----------
    container: A docker.models.container.Container object
    local_paths: Local paths of each file to copy

    Returns
    -------
    Returns a list of files in the docker container after the copying is done
    """

    if isinstance(local_paths, str):
        local_paths = [local_paths]

    for local_path in local_paths:
        docker_put(container, local_path)

    return (
        container.exec_run("ls", stderr=True, stdout=True).output.decode().split("\n")
    )


def docker_get_new_files(
    container: Container, output_dir: str | Path, ignore_files: list[str | Path]
):
    """
    Function to copy new files out of a container.
    All files in the work directory of 'container'
    will be copied out to 'output_dir',
    unless they appear in the 'ignore_files' list.

    The normal procedure is to run this in tandem with docker_batch_put(),
    so that only new files
    are copied out of the container. For example:

        ignore_files = docker_batch_put(
            container=container,
            local_paths=list_of_files_to_copy_into_container
        )

        container.exec_run(some_docker_command)

        docker_get_new_files(
            container=container,
            output_dir=output_dir,
            ignore_files=ignore_files
        )

    Parameters
    ----------
    container: A docker.models.container.Container object
    output_dir: A local directory to save the output files to.
    ignore_files: List of files to ignore (i.e to not copy)

    Returns
    -------
    """

    new_files = [
        x
        for x in container.exec_run("ls", stderr=True, stdout=True)
        .output.decode()
        .split("\n")
        if x not in ignore_files
    ]

    # Make output directory if it doesn't exist

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Collect output files

    for output_file in new_files:
        output_path = output_dir.joinpath(output_file)

        docker_get(container, output_path)

        if output_path.exists():
            logger.debug(f"Saved to {output_path}")
        else:
            raise FileNotFoundError(f"Unable to find {output_path}")
