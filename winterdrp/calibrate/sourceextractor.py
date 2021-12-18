import os.path
import docker
import logging
import io
import tarfile
from pathlib import Path
from winterdrp.paths import calibration_config_dir

logger = logging.getLogger(__name__)

from docker.errors import DockerException

docker_dir = "/root/"


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


docker_image_name = "robertdstein/astrodocker"


def new_container(pull=False):
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


def run_sextractor(
        image: str,
        output_dir: str,
        output_files: list = None,
        config: str = os.path.join(calibration_config_dir, 'astrom.sex'),
        param: str = os.path.join(calibration_config_dir, 'astrom.param'),
        filter_name: str = os.path.join(calibration_config_dir, 'default.conv'),
        star_nnw: str = os.path.join(calibration_config_dir, 'default.nnw'),
        weights: str = None
):

    image_name = Path(image).stem

    container = new_container()
    container.attach()

    container.start()

    # Add image

    for path in [image, config, param, filter_name, star_nnw]:
        docker_put(container, path)

    output_catalog = docker_path(f'{image_name}.cat')

    cmd = f"/usr/bin/source-extractor {docker_path(image)} " \
          f"-c {docker_path(config)} " \
          f"-CATALOG_NAME {output_catalog} " \
          f"-PARAMETERS_NAME {docker_path(param)} " \
          f"-FILTER_NAME {docker_path(filter_name)} " \
          f"-STARNNW_NAME {docker_path(star_nnw)} " \
          f"-VERBOSE_TYPE QUIET "

    if weights is None:
        cmd += "-WEIGHT_TYPE None"
    else:
        docker_put(container, weights)

    # Add config file

    # cmd += "-c sex.config"

    # docker_put(container, config)
    # cmd += f"-c {os.path.basename(config)}"

    # Run sextractor

    log = container.exec_run(cmd, stderr=True, stdout=True)

    if not log.output == b"":
        logger.warning(f"Sextractor warning: {log.output.decode()}")

    # Collect output files

    if output_files is None:
        output_files = [os.path.basename(output_catalog)]

    for output_file in output_files:

        output_path = os.path.join(output_dir, output_file)

        docker_get(container, output_path)

        if os.path.exists(output_path):
            logger.debug(f"Saved to {output_path}")
        else:
            raise FileNotFoundError(f"Unable to find {output_path}")

    container.kill()
    container.remove()


run_sextractor(
    "/Users/robertstein/Data/WIRC/20200929/redux/image0262.fits",
    output_dir="/Users/robertstein/Data/WIRC/20200929/sextractor/"
)