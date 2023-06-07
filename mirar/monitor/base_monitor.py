"""
Script containing the :class:`~mirar.monitor.base_monitor.Monitor` class,
used for processing data in real time.
"""
import copy
import logging
import os
import sys
import threading
import time
from pathlib import Path
from queue import Queue
from threading import Thread
from typing import Optional

import numpy as np
from astropy import units as u
from astropy.time import Time
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer

from mirar.data import Dataset, Image, ImageBatch
from mirar.errors import ErrorReport, ErrorStack, ImageNotFoundError
from mirar.io import check_file_is_complete
from mirar.paths import (
    MONITOR_EMAIL_KEY,
    MONITOR_RECIPIENT_KEY,
    PACKAGE_NAME,
    RAW_IMG_SUB_DIR,
    __version__,
    base_raw_dir,
    get_output_path,
    max_n_cpu,
    raw_img_dir,
)
from mirar.pipelines import get_pipeline
from mirar.processors.csvlog import CSVLog
from mirar.processors.utils.cal_hunter import (
    CalHunter,
    CalRequirement,
    find_required_cals,
    update_requirements,
)
from mirar.processors.utils.image_loader import ImageLoader
from mirar.utils.send_email import send_gmail

logger = logging.getLogger(__name__)


class NewImageHandler(FileSystemEventHandler):
    """Class to watch a directory, and add newly-created files to a queue."""

    def __init__(self, queue):
        FileSystemEventHandler.__init__(self)
        self.queue = queue

    def on_created(self, event):
        if event.event_type == "created":
            self.queue.put(event)


class Monitor:
    """Class to 'monitor' a directory, watching for newly created files.
    It then reduces these files. It will watch for a fixed duration,
    and run a postprocessing step at some configurable time after starting.
    It can send automated email notifications.
    """

    def __init__(
        self,
        night: str,
        pipeline: str,
        cal_requirements: Optional[list[CalRequirement]] = None,
        realtime_configurations: str | list[str] = "default",
        postprocess_configurations: Optional[str | list[str]] = None,
        email_sender: Optional[str] = os.getenv(MONITOR_EMAIL_KEY),
        email_recipients: Optional[str | list] = os.getenv(MONITOR_RECIPIENT_KEY),
        midway_postprocess_hours: float = 16.0,
        final_postprocess_hours: float = 48.0,
        log_level: str = "INFO",
        raw_dir: str = RAW_IMG_SUB_DIR,
        base_raw_img_dir: Path = base_raw_dir,
    ):
        logger.info(f"Software version: {PACKAGE_NAME}=={__version__}")

        self.errorstack = ErrorStack()
        self.night = night
        self.pipeline_name = pipeline

        if not isinstance(realtime_configurations, list):
            realtime_configurations = [realtime_configurations]
        self.realtime_configurations = realtime_configurations

        self.postprocess_configurations = postprocess_configurations

        self.pipeline = get_pipeline(
            pipeline, night=night, selected_configurations=realtime_configurations
        )

        self.raw_image_directory = Path(
            raw_img_dir(
                sub_dir=self.pipeline.night_sub_dir,
                img_sub_dir=raw_dir,
                raw_dir=base_raw_img_dir,
            )
        )

        self.raw_image_directory.mkdir(parents=True, exist_ok=True)

        self.sub_dir = raw_dir

        self.log_level = log_level
        self.log_path = self.configure_logs(log_level)
        self.error_path = self.pipeline.get_error_output_path()

        self.final_postprocess_hours = float(final_postprocess_hours) * u.hour
        logger.info(f"Will terminate after {final_postprocess_hours} hours.")
        self.t_start = Time.now()

        self.midway_postprocess_hours = float(midway_postprocess_hours) * u.hour

        if self.midway_postprocess_hours > self.final_postprocess_hours:
            logger.warning(
                f"Midway postprocessing was set to {self.midway_postprocess_hours}, "
                "but the monitor has a shorter termination period of "
                f"{self.final_postprocess_hours}. Setting to to 95% of max wait."
            )
            self.midway_postprocess_hours = 0.95 * self.final_postprocess_hours

        check_email = np.sum([x is not None for x in [email_recipients, email_sender]])
        if np.sum(check_email) == 1:
            err = (
                "In order to send emails, you must specify both a sender"
                f" and a recipient. \n In this case, sender is {email_sender} "
                f"and recipient is {email_recipients}."
            )
            logger.error(err)
            raise ValueError(err)

        if np.sum(check_email) == 2:
            logger.info(
                f"Will send an email summary after "
                f"{self.midway_postprocess_hours} hours."
            )
            self.email_info = (email_sender, email_recipients)
            self.email_to_send = True

        else:
            logger.info("No email notification configured.")
            self.email_info = None
            self.email_to_send = False

        self.midway_postprocess_complete = False
        self.latest_csv_log = None

        self.processed_science_images = []
        self.processed_cal_images = []
        self.failed_images = []

        # default to "pipeline default cal requirements"

        if cal_requirements is None:
            cal_requirements = self.pipeline.default_cal_requirements

        self.archival_cals = ImageBatch()
        self.new_cals = ImageBatch()
        self.cal_requirements = copy.deepcopy(cal_requirements)

        if cal_requirements is not None:
            try:
                self.archival_cals = find_required_cals(
                    latest_dir=str(self.raw_image_directory),
                    night=night,
                    open_f=self.pipeline.load_raw_image,
                    requirements=cal_requirements,
                )
            except ImageNotFoundError as exc:
                err = "No CalHunter images found. Will need to rely on nightly data."
                logger.error(err)
                self.errorstack.add_report(
                    ErrorReport(
                        error=exc, processor_name=CalHunter.__name__, contents=[]
                    )
                )

    def get_cals(self) -> ImageBatch:
        """
        Returns a copy of the calibration images (new and archival)

        :return:
        """
        return copy.deepcopy(self.new_cals + self.archival_cals)

    def update_cals(self, new_calibration_image: Image):
        """
        Updates the calibration images by adding a new calibration image.
        The archival cal images are then rechecked, and only those which are still
        required are loaded.

        :param new_calibration_image: new image
        :return: None
        """
        self.new_cals.append(new_calibration_image)
        cal_requirements = copy.deepcopy(self.cal_requirements)
        cal_requirements = [
            x
            for x in update_requirements(cal_requirements, self.new_cals)
            if not x.success
        ]

        cal_requirements = update_requirements(cal_requirements, self.archival_cals)
        new_archival_cals = ImageBatch()

        for archival_cal in self.archival_cals:
            for req in cal_requirements:
                for batch in req.data.values():
                    if archival_cal in batch:
                        if archival_cal not in new_archival_cals:
                            new_archival_cals.append(archival_cal)

        self.archival_cals = new_archival_cals

    def summarise_errors(
        self,
        errorstack: ErrorStack,
    ):
        """Create a text summary using an errorstack and the list
        of processed images. Sends an email of this if configured
        to do so, or prints otherwise.

        :param errorstack: list of errors to summarise
        :return: None
        """

        error_summary = errorstack.summarise_error_stack(verbose=False)
        summary = (
            f"Processed a total of {len(self.processed_science_images)}"
            f" science images. \n\n {error_summary} \n"
        )

        logger.info(f"Writing error log to {self.error_path}")
        errorstack.summarise_error_stack(verbose=True, output_path=self.error_path)

        if self.email_info is not None:
            sender, recipients = self.email_info

            subject = f"{self.pipeline_name}: Summary for night {self.night}"

            attachments = [self.log_path, self.error_path]

            # Send the latest CSV log if there is one
            if self.latest_csv_log is not None:
                attachments.append(self.latest_csv_log)

            send_gmail(
                email_sender=sender,
                email_recipients=recipients,
                email_subject=subject,
                email_text=summary,
                attachments=attachments,
            )
        else:
            print(summary)

    def configure_logs(self, log_level="INFO"):
        """Function to configure the log level for the python logger.
        Posts the log to the terminal and also writes it to a file.

        :param log_level: python log level
        :return: lof file path
        """

        log_output_path = get_output_path(
            base_name=f"{self.night}_processing_log.txt",
            dir_root=self.pipeline.night_sub_dir,
        )

        try:
            os.makedirs(os.path.dirname(log_output_path))
        except OSError:
            pass

        log = logging.getLogger("mirar")

        handler = logging.FileHandler(log_output_path)
        # handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter(
            "%(asctime)s: %(name)s [l %(lineno)d] - %(levelname)s - %(message)s"
        )
        handler.setFormatter(formatter)
        log.addHandler(handler)
        log.setLevel(log_level)

        root = logging.getLogger()
        root.setLevel(log_level)

        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(log_level)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        handler.setFormatter(formatter)
        root.addHandler(handler)

        logger.info(f"Logging level: {self.log_level}, saving log to {log_output_path}")
        return log_output_path

    def process_realtime(self):
        """Function to initiate the actual monitoring.

        :return: None
        """
        # create queue
        monitor_queue = Queue()

        workers = []

        n_cpu = max_n_cpu

        for _ in range(n_cpu):
            # Set up a worker thread to process database load
            worker = Thread(target=self.process_load_queue, args=(monitor_queue,))
            worker.daemon = True
            worker.start()

            workers.append(worker)

        # setup watchdog to monitor directory for trigger files
        logger.info(f"Watching {self.raw_image_directory}")

        event_handler = NewImageHandler(monitor_queue)
        observer = Observer()
        observer.schedule(event_handler, path=str(self.raw_image_directory))
        observer.start()

        try:
            while (Time.now() - self.t_start) < self.final_postprocess_hours:
                time.sleep(2)
        finally:
            logger.info("No longer waiting for new images.")
            observer.stop()
            observer.join()
            self.postprocess()

    def update_error_log(self):
        """Function to overwrite the error file with the latest version.
        The error summary is cumulative, so this just updates the file.
        """
        self.errorstack.summarise_error_stack(verbose=True, output_path=self.error_path)

    def postprocess(self):
        """Function to be run after some realtime postprocessing has been run.
        This function is called once after a configurable number of hours
        (typically when the data is expected to be done), and then again
        when the monitor stops watching the directory.

        :return: None
        """
        self.update_error_log()

        logger.info("Running postprocess steps")

        if self.postprocess_configurations is not None:
            postprocess_config = [
                ImageLoader(
                    load_image=self.pipeline.unpack_raw_image,
                    input_sub_dir=self.sub_dir,
                    input_img_dir=str(Path(self.raw_image_directory)).split(
                        self.pipeline_name, maxsplit=1
                    )[0],
                )
            ]

            postprocess_config += self.pipeline.postprocess_configuration(
                errorstack=self.errorstack,
                processed_images=[
                    os.path.basename(x) for x in self.processed_science_images
                ],
                selected_configurations=self.postprocess_configurations,
            )

            protected_key = "_monitor"
            while protected_key in self.pipeline.all_pipeline_configurations.keys():
                protected_key += "_2"

            self.pipeline.add_configuration(protected_key, postprocess_config)
            self.pipeline.set_configuration(protected_key)

            for processor in self.pipeline.all_pipeline_configurations[protected_key]:
                if isinstance(processor, CSVLog):
                    self.latest_csv_log = processor.get_output_path()

            _, errorstack = self.pipeline.reduce_images(
                dataset=Dataset(ImageBatch()),
                selected_configurations=protected_key,
                catch_all_errors=True,
            )
            self.errorstack += errorstack
            self.update_error_log()

    def process_load_queue(self, queue: Queue):
        """This is the worker thread function. It is run as a daemon
        threads that only exit when the main thread ends.

        Args
        ==========
          queue:  Queue() object
        """
        while True:
            if Time.now() - self.t_start > self.midway_postprocess_hours:
                if not self.midway_postprocess_complete:
                    self.midway_postprocess_complete = True
                    logger.info("Postprocess time!")
                    self.postprocess()
                    if self.email_to_send:
                        logger.info(
                            f"More than {self.midway_postprocess_hours} "
                            f"hours have elapsed. Sending summary email."
                        )
                        self.summarise_errors(errorstack=self.errorstack)

            if not queue.empty():
                event = queue.get()

                if event.src_path[-5:] == ".fits":
                    # Verify that file transfer is complete, useful for rsync latency

                    transfer_done = False

                    while not transfer_done:
                        transfer_done = check_file_is_complete(event.src_path)

                        if not transfer_done:
                            print(
                                "Seems like the file is not fully transferred. "
                                "Waiting a couple of seconds before trying again."
                            )
                            time.sleep(3)

                    try:
                        img = self.pipeline.load_raw_image(event.src_path)

                        is_science = img["OBSCLASS"] == "science"

                        if not is_science:
                            self.update_cals(img)

                        all_img = ImageBatch(img) + self.get_cals()

                        print(
                            f"Reducing {event.src_path} "
                            f"on thread {threading.get_ident()}, "
                            f"(science={is_science})"
                        )
                        _, errorstack = self.pipeline.reduce_images(
                            dataset=Dataset(all_img),
                            selected_configurations=self.realtime_configurations,
                            catch_all_errors=True,
                        )
                        self.errorstack += errorstack
                        self.update_error_log()

                        if is_science:
                            self.processed_science_images.append(event.src_path)
                        else:
                            self.processed_cal_images.append(event.src_path)

                    # RS: Please forgive me for this coding sin
                    # I just want the monitor to never crash
                    except Exception as exc:  # pylint: disable=broad-except
                        err_report = ErrorReport(
                            exc, "monitor", contents=[event.src_path]
                        )
                        self.errorstack.add_report(err_report)
                        self.update_error_log()
                        self.failed_images.append(event.src_path)

            else:
                time.sleep(1)
