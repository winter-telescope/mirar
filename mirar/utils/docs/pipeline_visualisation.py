"""
Module for generating visualisations of
:class:`~mirar.pipeline.base_pipeline.Pipeline` objects.
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt

from mirar.paths import doc_dir
from mirar.pipelines import Pipeline, get_pipeline
from mirar.processors.base_processor import (
    BaseImageProcessor,
    BaseProcessor,
    BaseSourceGenerator,
    BaseSourceProcessor,
)

logger = logging.getLogger(__name__)

docs_source_dir = doc_dir / "source"

base_autogen_dir = docs_source_dir / "autogen"
base_autogen_dir.mkdir(parents=True, exist_ok=True)

flowchart_dir = base_autogen_dir / "flowcharts"
flowchart_dir.mkdir(exist_ok=True)


def get_save_path(pipeline: str, configs: str) -> Path:
    """
    Get output save path for a pipeline diagram

    :param pipeline: Pipeline ised
    :param configs: Configs used
    :return: path to save
    """
    return flowchart_dir.joinpath(f"{pipeline}/{configs}.png")


def flowify(
    processor_list: list[BaseProcessor],
    output_path: Path,
    include_stats: bool = False,
):
    """
    Function to generate a diagram summarising all
    :class:`~wintedrp.processors.BaseProcessor` objects
    in a given pipeline configuration

    :param processor_list: list of processors to visualise
    :param output_path: Path to save diagram
    :param include_stats: Include statistics in diagram
    :return: None
    """

    plt.figure(
        figsize=(12.0 + 6.0 * include_stats, 2.0 + 0.35 * len(processor_list)),
        dpi=300.0,
    )
    plt.subplot(111)

    y_scale = 1.0 / float(len(processor_list))

    base_offset = 0.8
    x_offset_name = 0.1
    x_offset_description = 0.6

    x_offset_stats = 1.1
    x_offset_errors = 1.25

    for i, processor in enumerate(processor_list):
        y_0 = 1.0 - y_scale * (i + base_offset + 0.7)
        y_1 = 1.0 - y_scale * (i + base_offset)

        xy_name = (x_offset_name, y_0)
        xytext = (x_offset_name, y_1)
        xy_description = (x_offset_description, y_0)
        xytext_description = (x_offset_description, y_1)

        if isinstance(processor, BaseImageProcessor):
            class_kwargs = {"color": "g"}
            blocks = "images"
        elif isinstance(processor, BaseSourceGenerator):
            class_kwargs = {"color": "purple"}
            blocks = "images"
        elif isinstance(processor, BaseSourceProcessor):
            class_kwargs = {"color": "red"}
            blocks = "sourcetables"
        else:
            raise ValueError(f"processor type ({type(processor)} not recognised")

        if i < len(processor_list) - 1:
            arrowprops = {
                "arrowstyle": "simple",
            }
        else:
            arrowprops = None

        annotate_args = {
            "xycoords": "axes fraction",
            "ha": "center",
            "arrowprops": arrowprops,
            "bbox": {
                "boxstyle": "round",
                "facecolor": "wheat",
                "alpha": 0.5,
            },
        }

        plt.annotate(
            text=processor.__class__.__name__,
            xy=xy_name,
            xytext=xytext,
            **annotate_args,
            **class_kwargs,
        )

        plt.annotate(
            text=processor.description(),
            xy=xy_description,
            xytext=xytext_description,
            **annotate_args,
            **class_kwargs,
        )

        if include_stats:

            err_stack = processor.latest_error_stack

            msg = (
                f"{processor.latest_n_input_blocks} {blocks}, "
                f"{processor.latest_n_input_batches} batches, "
                f"{len(err_stack.reports)} errors"
            )

            plt.annotate(
                text=msg,
                xy=(x_offset_stats, y_0),
                xytext=(x_offset_stats, y_1),
                **annotate_args,
                **class_kwargs,
            )

            if len(err_stack.reports) > 0:
                err_names = [err.get_error_name() for err in err_stack.reports]
                err_counts = {name: err_names.count(name) for name in set(err_names)}
                err_msg = ", ".join([f"{v}x{k}" for k, v in err_counts.items()])

                annotate_args["ha"] = "left"
                annotate_args["arrowprops"] = {
                    "arrowstyle": "<-",
                }

                plt.annotate(
                    text=err_msg,
                    xy=(x_offset_stats + 0.1, y_1),
                    xytext=(x_offset_errors, y_1),
                    **annotate_args,
                    **class_kwargs,
                )

    logger.info(f"Saving to {output_path}")

    if not output_path.parent.exists():
        output_path.parent.mkdir(parents=True)

    plt.tight_layout()
    plt.axis("off")

    plt.savefig(output_path)
    plt.close()


def iterate_flowify(
    config: Optional[str | list[str]] = None,
    pipelines: Optional[str | list[str]] = None,
):
    """
    Function to iterate the visualisation of all configurations and pipelines

    :param config: config(s) to visualise (default of all)
    :param pipelines: pipeline(s) to visualise
    :return: None
    """
    if pipelines is None:
        pipelines = Pipeline.pipelines.keys()
    elif not isinstance(pipelines, list):
        pipelines = [pipelines]

    if config is None:
        selected_config = "default"
    else:
        selected_config = config

    for pipeline in pipelines:
        pipe = get_pipeline(
            pipeline,
            selected_configurations=selected_config,
            night=str(datetime.now()).split(" ", maxsplit=1)[0].replace("-", ""),
        )

        logger.info(f"Visualising {pipeline} pipeline")

        if config is None:
            config_list = pipe.all_pipeline_configurations.keys()
        else:
            config_list = pipe.selected_configurations

        for single_config in config_list:
            logger.info(f"Visualising {single_config} configuration")
            flowify(
                pipe.set_configuration(single_config),
                get_save_path(pipeline, single_config),
            )


if __name__ == "__main__":
    iterate_flowify()
