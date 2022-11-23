import os
from watchdog.events import FileSystemEventHandler
import threading
from threading import Thread
from queue import Queue
import time
import sys
from watchdog.observers import Observer
from winterdrp.pipelines import get_pipeline, PipelineConfigError
from winterdrp.errors import ErrorStack, ErrorReport
from winterdrp.utils.send_email import send_gmail
from winterdrp.paths import get_output_path, raw_img_dir, raw_img_sub_dir, __version__, package_name, base_raw_dir, \
    watchdog_email_key, watchdog_recipient_key
import numpy as np
import logging
from astropy.time import Time
from astropy import units as u
from winterdrp.processors.utils.cal_hunter import CalRequirement, find_required_cals
from pathlib import Path
import copy
from astropy.io import fits
from warnings import catch_warnings
import warnings
from astropy.utils.exceptions import AstropyUserWarning
from winterdrp.processors.utils.image_loader import ImageLoader
from winterdrp.processors.csvlog import CSVLog
from winterdrp.data import ImageBatch, Dataset

logger = logging.getLogger(__name__)


class NewImageHandler(FileSystemEventHandler):

    def __init__(self, queue):
        FileSystemEventHandler.__init__(self)
        self.queue = queue

    def on_created(self, event):
        if event.event_type == "created":
            self.queue.put(event)


class Monitor:

    def __init__(
            self,
            night: str,
            pipeline: str,
            cal_requirements: list[CalRequirement] = None,
            realtime_configurations: str | list[str] = "default",
            postprocess_configurations: str | list[str] = None,
            email_sender: str = os.getenv(watchdog_email_key),
            email_recipients: str | list = os.getenv(watchdog_recipient_key),
            midway_postprocess_hours: float = 16.,
            final_postprocess_hours: float = 48.,
            log_level: str = "INFO",
            raw_dir: str = raw_img_sub_dir,
            base_raw_img_dir: str = base_raw_dir
    ):

        logger.info(f"Software version: {package_name}=={__version__}")

        self.errorstack = ErrorStack()
        self.night = night
        self.pipeline_name = pipeline

        if not isinstance(realtime_configurations, list):
            realtime_configurations = [realtime_configurations]
        self.realtime_configurations = realtime_configurations

        self.postprocess_configurations = postprocess_configurations

        self.pipeline = get_pipeline(pipeline, night=night, selected_configurations=realtime_configurations)

        self.raw_image_directory = Path(raw_img_dir(
            sub_dir=self.pipeline.night_sub_dir,
            img_sub_dir=raw_dir,
            raw_dir=base_raw_img_dir
        ))

        self.sub_dir = raw_dir

        if not self.raw_image_directory.exists():
            for x in self.raw_image_directory.parents[::-1]:
                x.mkdir(exist_ok=True)
            self.raw_image_directory.mkdir()

        self.log_level = log_level
        self.log_path = self.configure_logs(log_level)
        self.error_path = self.get_error_output_path()

        check_email = np.sum([x is not None for x in [email_recipients, email_sender]])
        if np.sum(check_email) == 1:
            err = "In order to send emails, you must specify both a a sender and a recipient. \n" \
                  f"In this case, sender is {email_sender} and recipient is {email_recipients}."
            logger.error(err)
            raise ValueError(err)

        self.final_postprocess_hours = float(final_postprocess_hours) * u.hour
        logger.info(f"Will terminate after {final_postprocess_hours} hours.")
        self.t_start = Time.now()

        self.midway_postprocess_hours = float(midway_postprocess_hours) * u.hour

        if self.midway_postprocess_hours > self.final_postprocess_hours:
            logger.warning(f"Midway postprocessing was set to {self.midway_postprocess_hours}, "
                           f"but the monitor has a shorter termination period of {self.final_postprocess_hours}. "
                           f"Setting to to 95% of max wait.")
            self.midway_postprocess_hours = 0.95 * self.final_postprocess_hours

        if np.sum(check_email) == 2:
            logger.info(f"Will send an email summary after {self.midway_postprocess_hours} hours.")
            self.email_info = (email_sender, email_recipients)
            self.email_to_send = True

        else:
            logger.info("No email notification configured.")
            self.email_info = None
            self.email_to_send = False

        self.midway_postprocess_complete = False
        self.latest_csv_log = None

        self.processed_science_images = []
        self.corrupted_images = []

        # default to "pipeline default cal requirements"

        if cal_requirements is None:
            cal_requirements = self.pipeline.default_cal_requirements

        if cal_requirements is not None:
            self.cal_images = find_required_cals(
                latest_dir=str(self.raw_image_directory),
                night=night,
                open_f=self.pipeline.load_raw_image,
                requirements=cal_requirements
            )
        else:
            self.cal_images = ImageBatch()

    def summarise_errors(
            self,
            errorstack: ErrorStack,
    ):

        error_summary = errorstack.summarise_error_stack(verbose=False)
        summary = f"Processed a total of {len(self.processed_science_images)} science images. " \
                  f"\n\n" + error_summary + " \n"

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
                attachments=attachments
            )
        else:
            print(summary)

    def process_full_night(
            self,
    ):
        batches, errorstack = self.pipeline.reduce_images(Dataset(), catch_all_errors=True)
        self.summarise_errors(errorstack=errorstack)

    def get_error_output_path(self) -> str:
        error_output_path = get_output_path(
            base_name=f"{self.night}_error_stack.txt",
            dir_root=self.pipeline.night_sub_dir,
        )

        try:
            os.makedirs(os.path.dirname(error_output_path))
        except OSError:
            pass

        return error_output_path

    def configure_logs(self, log_level="INFO"):

        log_output_path = get_output_path(
            base_name=f"{self.night}_processing_log.txt",
            dir_root=self.pipeline.night_sub_dir,
        )

        try:
            os.makedirs(os.path.dirname(log_output_path))
        except OSError:
            pass

        log = logging.getLogger("winterdrp")

        handler = logging.FileHandler(log_output_path)
        # handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter('%(asctime)s: %(name)s [l %(lineno)d] - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        log.addHandler(handler)
        log.setLevel(log_level)

        root = logging.getLogger()
        root.setLevel(log_level)

        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(log_level)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        root.addHandler(handler)

        logger.info(f"Logging level: {self.log_level}, saving log to {log_output_path}")
        return log_output_path

    def process_realtime(self):
        # create queue
        watchdog_queue = Queue()

        workers = []

        n_cpu = max(1, int(os.cpu_count() / 2))

        for i in range(n_cpu):
            # Set up a worker thread to process database load
            worker = Thread(target=self.process_load_queue, args=(watchdog_queue,))
            worker.daemon = True
            worker.start()

            workers.append(worker)

        # setup watchdog to monitor directory for trigger files
        logger.info(f"Watching {self.raw_image_directory}")

        event_handler = NewImageHandler(watchdog_queue)
        observer = Observer()
        observer.schedule(event_handler, path=str(self.raw_image_directory))
        observer.start()

        try:
            while (Time.now() - self.t_start) < self.final_postprocess_hours:
                time.sleep(2)
        finally:
            logger.info(f"No longer waiting for new images.")
            observer.stop()
            observer.join()
            self.postprocess()

    def update_error_log(self):
        self.errorstack.summarise_error_stack(verbose=True, output_path=self.error_path)

    def postprocess(self):
        self.update_error_log()

        logger.info("Running postprocess steps")

        if self.postprocess_configurations is not None:

            postprocess_config = [ImageLoader(
                load_image=self.pipeline.unpack_raw_image,
                input_sub_dir=self.sub_dir,
                input_img_dir=str(Path(self.raw_image_directory)).split(self.pipeline_name)[0]
            )]

            postprocess_config += self.pipeline.postprocess_configuration(
                errorstack=self.errorstack,
                processed_images=[os.path.basename(x) for x in self.processed_science_images],
                selected_configurations=self.postprocess_configurations
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
                catch_all_errors=True
            )
            self.errorstack += errorstack
            self.update_error_log()

    def process_load_queue(self, q):
        '''This is the worker thread function. It is run as a daemon
           threads that only exit when the main thread ends.

           Args
           ==========
             q:  Queue() object
        '''
        while True:

            if Time.now() - self.t_start > self.midway_postprocess_hours:
                if not self.midway_postprocess_complete:
                    self.midway_postprocess_complete = True
                    logger.info("Postprocess time!")
                    self.postprocess()
                    if self.email_to_send:
                        logger.info(f"More than {self.midway_postprocess_hours} hours have elapsed. "
                                    f"Sending summary email.")
                        self.summarise_errors(errorstack=self.errorstack)

            if not q.empty():
                event = q.get()

                if event.src_path[-5:] == ".fits":

                    # Verify that file transfer is complete, useful for rsync latency

                    # Disclaimer: I (Robert) do not feel great about having written this code block
                    # It seems to works though, let's hope no one finds out
                    # I will cover my tracks by hiding the astropy warning which inspired this block,
                    # informing the user that the file is not as long as expected

                    check = False

                    while not check:
                        with catch_warnings():
                            warnings.filterwarnings('ignore', category=AstropyUserWarning)
                            try:
                                with fits.open(event.src_path) as hdul:
                                    check = hdul._file.tell() == hdul._file.size
                            except OSError:
                                pass

                            if not check:
                                print("Seems like the file is not fully transferred. "
                                      "Waiting a couple of seconds before trying again.")
                                time.sleep(3)

                    try:

                        img = self.pipeline.load_raw_image(event.src_path)

                        is_science = img["OBSCLASS"] == "science"

                        all_img = ImageBatch(img) + copy.deepcopy(self.cal_images)

                        if not is_science:
                            print(f"Skipping {event.src_path} (calibration image)")
                        else:
                            print(f"Reducing {event.src_path} (science image) on thread {threading.get_ident()}")
                            _, errorstack = self.pipeline.reduce_images(
                                dataset=Dataset(all_img),
                                selected_configurations=self.realtime_configurations,
                                catch_all_errors=True
                            )
                            self.processed_science_images.append(event.src_path)
                            self.errorstack += errorstack
                            self.update_error_log()

                    except Exception as e:
                        err_report = ErrorReport(e, "monitor", contents=[event.src_path])
                        self.errorstack.add_report(err_report)
                        self.update_error_log()

            else:
                time.sleep(1)
