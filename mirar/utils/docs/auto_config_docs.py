"""
Module to auto-generate rst files for each pipeline in the project.
"""

import logging
import textwrap
from datetime import datetime
from pathlib import Path

from rstcloth import RstCloth

from mirar.pipelines import Pipeline, get_pipeline
from mirar.utils.docs.pipeline_visualisation import (
    base_autogen_dir,
    get_save_path,
    iterate_flowify,
)

logger = logging.getLogger(__name__)


def get_rst_pipeline_path(pipeline: str) -> Path:
    """
    Get output save path for a pipeline diagram

    :param pipeline: Pipeline name
    :return: path to save
    """
    return base_autogen_dir.joinpath(f"{pipeline}.rst")


def get_rst_config_path(pipeline: str, config: str) -> Path:
    """
    Get output save path for a pipeline diagram

    :param pipeline: Pipeline ised
    :param config: Configs used
    :return: path to save
    """
    return base_autogen_dir.joinpath(f"{pipeline}.{config}.rst")


def auto_top_level_rst(pipelines: list[str]):
    """
    Function to generate a top-level rst file for all pipelines

    :return: None
    """

    output_path = base_autogen_dir.joinpath("pipelines.rst")

    with open(output_path, "w", encoding="utf8") as output_file:
        doc = RstCloth(output_file)
        doc.title("Pipeline Configurations")
        doc.newline()
        doc.h2("Pipelines")
        doc.content(
            textwrap.dedent(
                """This page contains a list of all available pipelines. 
            Each pipeline has a number of configurations, 
            and each configuration has a number of processors.

            Every configuration has an automatically-generated flowchart, 
            which can be found by following the links below.
            """
            )
        )
        doc.newline()

        doc.content(".. toctree::")
        doc.content("\t:maxdepth: 1")
        doc.content("\t:caption: Pipelines")
        doc.newline()

        for pipeline in sorted(pipelines):
            doc.content(f"\t{get_rst_pipeline_path(pipeline).stem}")

        doc.newline()


def auto_rst_pipeline(pipeline: str, configs: list[str]):
    """
    Function to generate a diagram summarising all
    :class:`~wintedrp.processors.BaseProcessor` objects
    in a given pipeline configuration

    :param pipeline: Pipeline name
    :param configs: List of configurations
    :return: None
    """

    output_path = get_rst_pipeline_path(pipeline)

    with open(output_path, "w", encoding="utf8") as output_file:
        doc = RstCloth(output_file)
        doc.title(f"{pipeline}")
        doc.newline()
        doc.content(
            textwrap.dedent(
                f"""This page contains a list of all {len(configs)}
            available configurations for {pipeline}. 

            Every configuration has an automatically-generated flowchart, 
            which can be found by clicking on the links below.
            """
            )
        )
        doc.newline()

        doc.content(".. toctree::")
        doc.content("    :maxdepth: 1")
        doc.content("    :caption: Configurations")
        doc.newline()

        for config in sorted(configs):
            doc.content(f"    {get_rst_config_path(pipeline, config).stem}")

        doc.newline()


def auto_rst_config(pipeline: str, config: str):
    """
    Function to generate a diagram summarising all
    :class:`~wintedrp.processors.BaseProcessor` objects
    in a given pipeline configuration

    :param pipeline: Pipeline name
    :param config: Configuration name
    :return: None
    """

    output_path = get_rst_config_path(pipeline, config)

    image_path = get_save_path(pipeline, config)

    if not image_path.exists():
        iterate_flowify(config=config, pipelines=[pipeline])

    relative_image_path = image_path.relative_to(base_autogen_dir)

    with open(output_path, "w", encoding="utf8") as output_file:
        doc = RstCloth(output_file)
        doc.title(f"{pipeline} - {config}")
        doc.newline()
        doc.h2("Flowchart")
        doc.content(
            textwrap.dedent(
                f"""This is the auto-generated flowchart for the '{config}' 
            config of the {pipeline} pipeline. 
            """
            )
        )
        doc.newline()
        doc.content(f".. image:: {relative_image_path}")
        doc.newline()


def iterate_rst_generation():
    """
    Function to iterate the creation of a separate rst file for each pipeline

    :return: None
    """
    pipelines = Pipeline.pipelines.keys()

    auto_top_level_rst(pipelines)

    for pipeline in pipelines:
        pipe = get_pipeline(
            pipeline,
            selected_configurations="default",
            night=str(datetime.now()).split(" ", maxsplit=1)[0].replace("-", ""),
        )

        config_list = pipe.all_pipeline_configurations.keys()

        auto_rst_pipeline(pipeline, config_list)

        for single_config in config_list:
            auto_rst_config(pipeline, single_config)


if __name__ == "__main__":
    iterate_rst_generation()
