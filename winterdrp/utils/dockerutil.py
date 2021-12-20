import os
import io
import tarfile
import docker
import logging
from docker.errors import DockerException
from docker.models.containers import Container

logger = logging.getLogger(__name__)

docker_image_name = "robertdstein/astrodocker"
docker_dir = "/usr/src/astrodocker"


def new_container():
    try:
        client = docker.from_env()
    except DockerException:
        err = "Unable to connect to Docker daemon. Have you installed Docker, and started a daemon? " \
              "Find out more at https://www.docker.com "
        logger.error(err)
        raise ConnectionError(err)

    if len(client.images.list(docker_image_name)) < 1:
        logger.info(f"Pulling docker image {docker_image_name}")
        client.images.pull(docker_image_name)

    return client.containers.run(docker_image_name, tty=True, detach=True)


def docker_path(file):
    return os.path.join(docker_dir, os.path.basename(file))


def docker_get(container, local_path):

    container_path = docker_path(local_path)

    with open(local_path, 'wb') as f:
        bits, stat = container.get_archive(container_path)
        for chunk in bits:
            f.write(chunk)


def docker_put(container, local_path):
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode='w|') as tar, open(local_path, 'rb') as f:
        info = tar.gettarinfo(fileobj=f)
        info.name = os.path.basename(local_path)
        tar.addfile(info, f)

    return container.put_archive(docker_dir, stream.getvalue())


def docker_batch_put(
        container: Container,
        local_paths: str | list = None
):
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

    return container.exec_run("ls", stderr=True, stdout=True).output.decode().split("\n")


def docker_get_new_files(
        container: Container,
        output_dir: str,
        ignore_files: list
):
    """

    Parameters
    ----------
    container
    output_dir: A local directory to save the output files to.
    ignore_files: List of files to ignore (i.e to not copy)

    Returns
    -------

    """

    new_files = [
        x for x in container.exec_run("ls", stderr=True, stdout=True).output.decode().split("\n")
        if x not in ignore_files
    ]

    # Make output directory

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    # Collect output files

    for output_file in new_files:

        output_path = os.path.join(output_dir, output_file)

        docker_get(container, output_path)

        if os.path.exists(output_path):
            logger.debug(f"Saved to {output_path}")
        else:
            raise FileNotFoundError(f"Unable to find {output_path}")
