import os
from winterdrp.paths import raw_img_dir


def download_via_ssh(
        server: str,
        base_dir: str,
        night: str | int,
        pipeline: str
):
    username = input(f"Please enter your username for {server}: \n")

    source_dir = f"{username}@{server}:{os.path.join(base_dir, night)}/"

    output_dir = raw_img_dir(os.path.join(pipeline, night))

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    cmd = f"rsync -a -v --include '*.fits' --exclude '*' {source_dir} {output_dir}"

    os.system(cmd)

