import logging
from winterdrp.utils import execute, ExecutionError

logger = logging.getLogger(__name__)

local_swarp = True


class SwarpExecutionError(ExecutionError):
    pass


def run_swarp(
        output_dir: str = ".",
        keyword_string: str = "",
        run_local: bool = local_swarp
):

    cmd = f"swarp {keyword_string}"

    if "RESAMPLE_DIR" not in keyword_string:
        cmd += f" -RESAMPLE_DIR {output_dir}"

    logger.debug(f"Executing '{cmd}'")

    try:
        execute(cmd, output_dir=output_dir, local=run_local)
    except ExecutionError as err:
        raise SwarpExecutionError(err)
