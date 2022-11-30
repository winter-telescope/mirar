import argparse
import logging
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from winterdrp.paths import doc_dir
from winterdrp.pipelines import get_pipeline, Pipeline
from winterdrp.processors.base_chain import BaseChain
from winterdrp.processors.base_processor import BaseImageProcessor, BaseCandidateGenerator, BaseDataframeProcessor

logger = logging.getLogger(__name__)


def get_save_path(pipeline: str, configs: str) -> Path:
    return doc_dir.joinpath(f"flowcharts/{pipeline}/{configs}.png")


def flowify(chain: BaseChain, output_path: Path):
    processor_list = chain.get_child_processors()
    plt.figure(figsize=(12., 2. + 0.3*len(processor_list)), dpi=300.)
    ax = plt.subplot(111)

    y_scale = 1./float(len(processor_list))

    base_offset = 0.8
    x_offset_name = 0.1
    x_offset_description = 0.6

    for i, processor in enumerate(processor_list):

        y_0 = 1. - y_scale*(i+base_offset + 0.7)
        y_1 = 1. - y_scale * (i + base_offset)

        xy_name = (x_offset_name, y_0)
        xytext = (x_offset_name, y_1)
        xy_description = (x_offset_description, y_0)
        xytext_description = (x_offset_description, y_1)

        if isinstance(processor, BaseImageProcessor):
            class_kwargs = {"color": "g"}
        elif isinstance(processor, BaseCandidateGenerator):
            class_kwargs = {"color": "purple"}
        elif isinstance(processor, BaseDataframeProcessor):
            class_kwargs = {"color": "red"}
        else:
            raise Exception(f"processor type ({type(processor)} not recognised")

        if i < len(processor_list) - 1:
            arrowprops = {
                "arrowstyle": "simple",
            }
        else:
            arrowprops = None

        annotate_args = {
            "xycoords": "axes fraction",
            "ha": 'center',
            "arrowprops": arrowprops,
            "bbox": {
                "boxstyle": 'round',
                "facecolor": 'wheat',
                "alpha": 0.5,
            }
        }

        plt.annotate(
            text=processor.__class__.__name__,
            xy=xy_name,
            xytext=xytext,
            **annotate_args,
            **class_kwargs
        )

        plt.annotate(
            text=str(processor),
            xy=xy_description,
            xytext=xytext_description,
            **annotate_args,
            **class_kwargs
        )

    logger.info(f"Saving to {output_path}")

    if not output_path.parent.exists():
        output_path.parent.mkdir(parents=True)

    plt.tight_layout()
    plt.axis('off')

    plt.savefig(output_path)


def iterate_flowify(
        config: str | list[str] = None,
        pipelines: str | list[str] = None
):
    if pipelines is None:
        pipelines = Pipeline.pipelines.keys()
    elif not isinstance(pipelines, list):
        pipelines = [args.pipeline]

    if config is None:
        gc = "default"
    else:
        gc = config

    for pipeline in pipelines:

        pipe = get_pipeline(
            pipeline,
            selected_configurations=gc,
            night=str(datetime.now()).split(" ")[0].replace("-", ""),
        )

        if config is None:
            config_list = pipe.all_pipeline_configurations.keys()
        else:
            config_list = pipe.selected_configurations

        for c in config_list:
            flowify(pipe.set_configuration(c), get_save_path(pipeline, c))


if __name__ == "__main__":

    logger.setLevel("INFO")

    parser = argparse.ArgumentParser(
        description="winterdrp: An automated image reduction pipeline, developed for WINTER"
    )

    parser.add_argument(
        "-p",
        "--pipeline",
        default=None,
        help="Pipeline to be used"
    )
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        help="Pipeline configuration to be used"
    )

    args = parser.parse_args()

    iterate_flowify(config=args.config, pipelines=args.pipeline)


