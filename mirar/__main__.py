"""
Main executable for mirar. You can execute the code from the terminal like:

.. codeblock:: bash
    python -m mirar -args...
"""

import argparse
import logging
import sys
import tempfile

from astropy import units as u
from astropy.time import Time

from mirar.data import cache
from mirar.monitor.base_monitor import Monitor
from mirar.paths import PACKAGE_NAME, RAW_IMG_SUB_DIR, TEMP_DIR
from mirar.pipelines import Pipeline, get_pipeline
from mirar.processors.utils import ImageLoader
from mirar.utils.docs.pipeline_visualisation import flowify

logger = logging.getLogger(__name__)


def main():
    """
    Main executable for mirar. You can execute the code from the terminal like:

    .. codeblock:: bash
        python -m mirar -args...
    """
    parser = argparse.ArgumentParser(
        description=f"{PACKAGE_NAME}: Modular Image Reduction and Analysis Resource"
    )

    parser.add_argument(
        "-n", "--night", default=None, help="Sub-directory to use in the data directory"
    )
    parser.add_argument("-p", "--pipeline", help="Pipeline to be used", required=True)
    parser.add_argument(
        "-c", "--config", default=None, help="Pipeline configuration to be used"
    )
    parser.add_argument(
        "-pc",
        "--postprocessconfig",
        default=None,
        help="Pipeline configuration to be used",
    )
    parser.add_argument(
        "--logfile",
        default=None,
        help="If a path is passed, all logs will be written to this file.",
    )
    parser.add_argument("--level", default="INFO", help="Python logging level")
    parser.add_argument(
        "-b", "--batch", default=None, help="Only process a specific image batch"
    )
    parser.add_argument(
        "--download",
        help="Download images from server",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--failfast", help="Fail on first error", action="store_true", default=False
    )

    parser.add_argument("-m", "--monitor", action="store_true", default=False)
    parser.add_argument(
        "--emailrecipients",
        default=None,
        help="Spaceless comma-separated values of email recipients",
    )
    parser.add_argument(
        "--emailsender",
        default=None,
        help="One email sender",
    )
    parser.add_argument(
        "--midwaypostprocesshours",
        default=16.0,
        help="Time, in hours, to wait before sending a summary email",
    )
    parser.add_argument(
        "--finalpostprocesshours",
        default=48.0,
        help="Time, in hours, to wait before ceasing monitoring for new images",
    )
    parser.add_argument(
        "--rawdir",
        default=RAW_IMG_SUB_DIR,
        help="Subdirectory to look in for raw images of a given night",
    )

    args = parser.parse_args()

    if args.download:
        Pipeline.pipelines[args.pipeline.lower()].download_raw_images_for_night(
            night=args.night
        )

        logger.info("Download complete")

    night = args.night

    with tempfile.TemporaryDirectory(dir=TEMP_DIR) as temp_dir_path:
        print(f"Using cache {temp_dir_path}")

        cache.set_cache_dir(temp_dir_path)

        if args.monitor:
            if args.emailrecipients is not None:
                email_recipients = args.emailrecipients.split(",")
            else:
                email_recipients = None

            config = args.config
            if config is None:
                config = "realtime"

            if night is None:
                ln = Time.now()
                night = str(ln).split(" ", maxsplit=1)[0].replace("-", "")

            monitor = Monitor(
                pipeline=args.pipeline,
                night=night,
                realtime_configurations=config,
                postprocess_configurations=(
                    args.postprocessconfig.split(",")
                    if args.postprocessconfig is not None
                    else None
                ),
                log_level=args.level,
                final_postprocess_hours=args.finalpostprocesshours,
                midway_postprocess_hours=args.midwaypostprocesshours,
                email_sender=args.emailsender,
                email_recipients=email_recipients,
                raw_dir=args.rawdir,
            )
            monitor.process_realtime()

        else:
            # Set up logging

            log = logging.getLogger("mirar")

            if args.logfile is None:
                handler = logging.StreamHandler(sys.stdout)
            else:
                handler = logging.FileHandler(args.logfile)

            formatter = logging.Formatter(
                "%(name)s [l %(lineno)d] - %(levelname)s - %(message)s"
            )
            handler.setFormatter(formatter)
            log.addHandler(handler)
            log.setLevel(args.level)

            config = args.config
            if config is None:
                config = "default"

            if night is None:
                ln = Time.now() - 1.0 * u.day
                night = str(ln).split(" ", maxsplit=1)[0].replace("-", "")

            pipe = get_pipeline(
                args.pipeline,
                selected_configurations=config,
                night=night,
            )

            _, errorstack = pipe.reduce_images(
                catch_all_errors=not args.failfast,
            )

            processors = pipe.get_latest_configuration()
            flowify(processors, pipe.get_flowchart_output_path(), include_stats=True)

            if args.postprocessconfig is not None:
                post_config = [
                    x
                    for x in pipe.set_configuration(config)
                    if isinstance(x, ImageLoader)
                ][:1]
                post_config += pipe.postprocess_configuration(
                    errorstack=errorstack,
                    selected_configurations=args.postprocessconfig.split(","),
                )

                protected_key = "_new_postprocess"

                pipe.add_configuration(protected_key, post_config)
                pipe.set_configuration(protected_key)

                _, new_errorstack, _ = pipe.reduce_images(
                    selected_configurations=protected_key,
                    catch_all_errors=True,
                )
                errorstack += new_errorstack

            print(errorstack.summarise_error_stack(verbose=False))

            errorstack.summarise_error_stack(
                output_path=pipe.get_error_output_path(), verbose=True
            )

            logger.info("End of mirar execution")


if __name__ == "__main__":
    main()
