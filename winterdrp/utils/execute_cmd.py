import logging
import subprocess
import os
import docker
import shutil
from docker.errors import DockerException
from winterdrp.utils.dockerutil import new_container, docker_path, docker_batch_put, docker_get_new_files

logger = logging.getLogger(__name__)


class ExecutionError(Exception):
    pass


def run_local(
        cmd: str,
        output_dir: str = "."
):
    """
    Function to run on local machine using subprocess, with error handling.

    After the specified 'cmd' command has been run, any newly-generated files
    will be copied out of the current directory to 'output_dir'

    Parameters
    ----------
    cmd: A string containing the command you want to use to run sextractor. An example would be:
        cmd = '/usr/bin/source-extractor image0001.fits -c sex.config'
    output_dir: A local directory to save the output files to.

    Returns
    -------

    """

    try:

        # See what files are in the directory beforehand

        ignore_files = subprocess.run("ls", check=True, capture_output=True).stdout.decode().split("\n")

        logger.debug(f"Ignoring files {ignore_files}")

        # Run sextractor

        rval = subprocess.run(cmd, check=True, capture_output=True, shell=True)

        logger.debug(f'Sextractor ran successfully on image {cmd.split(" ")[1]}')
        logger.debug(rval.stdout.decode())

        try:
            os.makedirs(output_dir)
        except OSError:
            pass

        # Move new files to output dir

        new_files = [
            x for x in subprocess.run("ls", check=True, capture_output=True).stdout.decode().split("\n")
            if x not in ignore_files
        ]

        current_dir = subprocess.run("pwd", check=True, capture_output=True).stdout.decode().strip()

        logger.info(f"The following new files were created: {new_files}")

        for file in new_files:

            current_path = os.path.join(current_dir, file)
            output_path = os.path.join(output_dir, file)

            logger.debug(f"Moving {current_path} to {output_path}")

            shutil.move(current_path, output_path)

    except subprocess.CalledProcessError as err:
        msg = f"Error found when running with command: \n \n '{err.cmd}' \n \n" \
              f"This yielded a return code of {err.returncode}. The following traceback was found: \n {err.stderr.decode()}"
        logger.error(msg)
        raise ExecutionError(msg)


def temp_config(
        config_path: str,
        output_dir: str
) -> str:
    basename = f"temp_{os.path.basename(config_path)}"
    return os.path.join(output_dir, basename)


def run_docker(
        cmd: str,
        output_dir: str = "."
):
    """Function to run a command via Docker. A container will be generated automatically,
    but a Docker server must be running first. You can start one via the Desktop application,
    or on the command line with `docker start'.

    After the specified 'cmd' command has been run, any newly-generated files
     will be copied out of the container to 'output_dir'

    Parameters
    ----------
    cmd: A string containing the base arguments you want to use to run sextractor. An example would be:
        cmd = 'image01.fits -c sex.config'
    output_dir: A local directory to save the output files to.

    Returns
    -------

    """

    container = new_container()

    try:

        container.attach()

        container.start()

        split = cmd.split(" -")

        # Reorganise the commands so that each '-x' argument is grouped together
        # Basically still work even if someone puts the filename in a weird place

        sorted_split = []

        for i, arg in enumerate(split):
            sep = arg.split(" ")
            sorted_split.append(" ".join(sep[:2]))
            if len(sep) > 2:
                sorted_split[0] += " " + " ".join(sep[2:])

        new_split = []

        # Loop over sextractor command, and
        # copy everything that looks like a file into container

        copy_list = []

        config_file = None

        for i, arg in enumerate(sorted_split):
            sep = arg.split(" ")

            if sep[0] == "c":
                config_file = sep[1]

            new = list(sep)

            for j, x in enumerate(sep):
                if os.path.isfile(x):
                    new[j] = docker_path(sep[j])
                    copy_list.append(sep[j])

            new_split.append(" ".join(new))

        cmd = " -".join(new_split)

        # Be extra clever: go through config file and check there too!

        new_config_file = []

        if config_file is not None:
            with open(config_file, "rb") as f:
                for line in f.readlines():
                    args = [x for x in line.decode().split(" ") if x not in [""]]
                    new_args = list(args)
                    for i, arg in enumerate(args):
                        if os.path.isfile(arg):
                            copy_list.append(arg)
                            new_args[i] = docker_path(arg)
                    new_config_file.append(" ".join(new_args))

            temp_config_path = temp_config(config_file, output_dir)

            with open(temp_config_path, "w") as g:
                g.writelines(new_config_file)

            copy_list.append(temp_config_path)

            cmd = cmd.replace(docker_path(config_file), docker_path(temp_config_path))

        # Copy in files, and see what files are already there

        copy_list = list(set(copy_list))

        logger.debug(f"Copying {copy_list} into container")

        ignore_files = docker_batch_put(
            container=container,
            local_paths=copy_list
        )

        # Run sextractor

        log = container.exec_run(cmd, stderr=True, stdout=True)

        if not log.output == b"":
            logger.info(f"Output: {log.output.decode()}")

        if not log.exit_code == 0:
            err = f"Error running command: \n '{cmd}'\n which resulted in returncode '{log.exit_code}' and" \
                  f"the following error message: \n '{log.output.decode()}'"
            logger.error(err)
            raise subprocess.CalledProcessError(
                returncode=log.exit_code,
                cmd=cmd,
                stderr=log.output.decode()
            )

        # Copy out any files which did not exist before running sextractor

        docker_get_new_files(
            container=container,
            output_dir=output_dir,
            ignore_files=ignore_files
        )

    except docker.errors.APIError as err:
        logger.error(err)
        raise ExecutionError(err)
    finally:
        # In any case, clean up by killing the container and removing files
        container.kill()
        container.remove()


def execute(
        cmd,
        output_dir=".",
        local=True
):
    logger.debug(f"Using '{['local', 'docker'][local]}' "
                 f" installation to run `{cmd}`")
    if local:
        run_local(cmd, output_dir=output_dir)
    else:
        run_docker(cmd, output_dir=output_dir)