import os.path
import docker
import logging
import io
import tarfile
from pathlib import Path

logger = logging.getLogger(__name__)


def docker_get(container, container_path, local_path):
    with open(local_path, 'wb') as f:
        bits, stat = container.get_archive(container_path)
        for chunk in bits:
            f.write(chunk)


def docker_put(container, local_path, container_dir):
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode='w|') as tar, open(local_path, 'rb') as f:
        info = tar.gettarinfo(fileobj=f)
        info.name = os.path.basename(local_path)
        tar.addfile(info, f)

    return container.put_archive(container_dir, stream.getvalue())


def new_container():
    client = docker.from_env()
    return client.containers.run('robertdstein/astrodocker:latest', tty=True, detach=True)


def run_sextractor(image, output_path, ):

    image_name = Path(image).stem

    container = new_container()
    container.attach()

    container_img_path = f"/root/{os.path.basename(image)}"

    container.start()

    docker_put(container, image, os.path.dirname(container_img_path))

    cmd = f"/usr/bin/source-extractor {container_img_path} -c sex.config -CATALOG_NAME {image_name}.cat"

    log = container.exec_run(cmd, stderr=True, stdout=True)
    if not log.output == b"":
        logger.error(log.output.decode())
        raise Exception(f"Runnning Sextractor should give no output,"
                        f"but instead returned the following output:{log.output.decode()}")

    docker_get(container, container_img_path, output_path)

    if os.path.exists(output_path):
        logger.debug(f"Saved to {output_path}")
    else:
        raise FileNotFoundError(f"Unable to find {output_path}")

    container.kill()
    container.remove()


run_sextractor(
    "/Users/robertstein/Data/WIRC/20200929/redux/image0266.fits",
    "/Users/robertstein/Data/WIRC/20200929/sextractor/image0266.cat"
)