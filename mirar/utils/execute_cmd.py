"""
Module for executing bash commands
"""
import logging
import os
import shutil
import subprocess
from pathlib import Path

import docker

from mirar.utils.dockerutil import (
    docker_batch_put,
    docker_get_new_files,
    docker_path,
    new_container,
)

logger = logging.getLogger(__name__)


class ExecutionError(Exception):
    """Error relating to executing bash command"""


DEFAULT_TIMEOUT = 300.0


def run_local(cmd: str, output_dir: str = ".", timeout: float = DEFAULT_TIMEOUT):
    """
    Function to run on local machine using subprocess, with error handling.

    After the specified 'cmd' command has been run, any newly-generated files
    will be copied out of the current directory to 'output_dir'

    Parameters
    ----------
    cmd: A string containing the command you want to use to run sextractor.
    An example would be:
        cmd = '/usr/bin/source-extractor image0001.fits -c sex.config'
    output_dir: A local directory to save the output files to.
    timeout: Time to timeout in seconds

    Returns
    -------

    """

    try:
        # See what files are in the directory beforehand

        ignore_files = (
            subprocess.run("ls", check=True, capture_output=True)
            .stdout.decode()
            .split("\n")
        )

        # Run sextractor

        rval = subprocess.run(
            cmd, check=True, capture_output=True, shell=True, timeout=timeout
        )

        msg = "Successfully executed command. "

        if rval.stdout.decode() != "":
            msg += f"Found the following output: {rval.stdout.decode()}"
        logger.debug(msg)

        try:
            os.makedirs(output_dir)
        except OSError:
            pass

        # Move new files to output dir

        new_files = [
            x
            for x in subprocess.run("ls", check=True, capture_output=True)
            .stdout.decode()
            .split("\n")
            if x not in ignore_files
        ]

        current_dir = (
            subprocess.run("pwd", check=True, capture_output=True)
            .stdout.decode()
            .strip()
        )

        if len(new_files) > 0:
            logger.debug(
                f"The following new files were created in the current directory: "
                f"{new_files}"
            )

        for file in new_files:
            current_path = os.path.join(current_dir, file)
            output_path = os.path.join(output_dir, file)

            logger.info(f"File saved to {output_path}")

            shutil.move(current_path, output_path)

    except subprocess.CalledProcessError as err:
        msg = (
            f"Error found when running with command: \n \n '{err.cmd}' \n \n"
            f"This yielded a return code of {err.returncode}. "
            f"The following traceback was found: \n {err.stderr.decode()}"
        )
        logger.error(msg)
        raise ExecutionError(msg) from err


def temp_config(config_path: str | Path, output_dir: str | Path) -> Path:
    """
    Get a

    :param config_path:
    :param output_dir:
    :return:
    """
    basename = f"temp_{Path(config_path).name}"
    return Path(output_dir).joinpath(basename)


def run_docker(cmd: str, output_dir: Path | str = "."):
    """Function to run a command via Docker.
    A container will be generated automatically,
    but a Docker server must be running first.
    You can start one via the Desktop application,
    or on the command line with `docker start'.

    After the specified 'cmd' command has been run, any newly-generated files
     will be copied out of the container to 'output_dir'

    Parameters
    ----------
    cmd: A string containing the base arguments you want to use to run sextractor.
    An example would be:
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
        # Go through everything that looks like a file with paths in it after

        copy_list = []
        temp_files = []

        files_of_files = []

        for i, arg in enumerate(sorted_split):
            sep = arg.split(" ")

            if sep[0] == "c":
                files_of_files.append(sep[1])

            new = list(sep)

            for j, x in enumerate(sep):
                if len(x) > 0:
                    if os.path.isfile(x):
                        new[j] = docker_path(sep[j])
                        copy_list.append(sep[j])
                    elif x[0] == "@":
                        files_of_files.append(x[1:])
                    elif os.path.isdir(os.path.dirname(x)):
                        new[j] = docker_path(sep[j])

            new_split.append(" ".join(new))

        cmd = " -".join(new_split)

        # Be extra clever: go through files and check there too!

        logger.debug(
            f"Found the following files which should contain paths: {files_of_files}"
        )

        for path in files_of_files:
            new_file = []

            with open(path, "rb", encoding="utf8") as local_file:
                for line in local_file.readlines():
                    args = [x for x in line.decode().split(" ") if x not in [""]]
                    new_args = list(args)
                    for i, arg in enumerate(args):
                        if os.path.isfile(arg):
                            copy_list.append(arg)
                            new_args[i] = docker_path(arg)
                        elif os.path.isfile(arg.strip("\n")):
                            copy_list.append(arg.strip("\n"))
                            new_args[i] = str(docker_path(arg.strip("\n"))) + "\n"
                    new_file.append(" ".join(new_args))

            temp_file_path = temp_config(path, output_dir)

            with open(temp_file_path, "w", encoding="utf8") as temp_file:
                temp_file.writelines(new_file)

            copy_list.append(temp_file_path)

            cmd = cmd.replace(path + " ", str(docker_path(temp_file_path)) + " ")

        # Copy in files, and see what files are already there

        copy_list = list(set(copy_list))

        logger.debug(f"Copying {copy_list} into container")

        ignore_files = docker_batch_put(container=container, local_paths=copy_list)

        # Run command

        log = container.exec_run(cmd, stderr=True, stdout=True)

        for temp_file_path in temp_files:
            logger.debug(f"Deleting temporary file {temp_file_path}")
            os.remove(temp_file_path)

        if not log.output == b"":
            logger.info(f"Output: {log.output.decode()}")

        if not log.exit_code == 0:
            err = (
                f"Error running command: \n '{cmd}'\n "
                f"which resulted in returncode '{log.exit_code}' and"
                f"the following error message: \n '{log.output.decode()}'"
            )
            logger.error(err)
            raise subprocess.CalledProcessError(
                returncode=log.exit_code, cmd=cmd, stderr=log.output.decode()
            )

        # Copy out any files which did not exist before running sextractor

        docker_get_new_files(
            container=container, output_dir=output_dir, ignore_files=ignore_files
        )

    except docker.errors.APIError as err:
        logger.error(err)
        raise ExecutionError(err) from err
    finally:
        # In any case, clean up by killing the container and removing files
        container.kill()
        container.remove()


def execute(
    cmd: str,
    output_dir: Path | str = ".",
    local: bool = True,
    timeout: float = DEFAULT_TIMEOUT,
):
    """
    Generically execute a command either via bash or a docker container

    :param cmd: command
    :param output_dir: output directory for command
    :param local: boolean whether use local or docker
    :param timeout: timeout for local execution
    :return: None
    """
    logger.debug(
        f"Using '{['docker', 'local'][local]}' " f" installation to run `{cmd}`"
    )
    if local:
        run_local(cmd, output_dir=output_dir, timeout=timeout)
    else:
        run_docker(cmd, output_dir=output_dir)
