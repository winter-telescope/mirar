import os
from watchdog.events import FileSystemEventHandler
from threading import Thread
from queue import Queue
import time
import sys
from watchdog.observers import Observer
from winterdrp.pipelines import get_pipeline
from winterdrp.errors import ErrorStack
from winterdrp.utils.send_email import send_gmail
from winterdrp.paths import get_output_path, raw_img_dir, raw_img_sub_dir
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
            realtime_configurations: str | list[str] = None,
            email_sender: str = None,
            email_recipients: str | list = None,
            email_wait_hours: float = 24.,
            max_wait_hours: float = 48.,
            log_level="INFO",
            raw_dir: str = raw_img_sub_dir
    ):

        self.errorstack = ErrorStack()
        self.night = night
        self.pipeline_name = pipeline

        if not isinstance(realtime_configurations, list):
            realtime_configurations = [realtime_configurations]
        self.realtime_configurations = realtime_configurations

        self.pipeline = get_pipeline(pipeline, night=night, selected_configurations=realtime_configurations)

        self.raw_image_directory = Path(raw_img_dir(sub_dir=self.pipeline.night_sub_dir, img_sub_dir=raw_dir))

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

        self.max_wait_hours = float(max_wait_hours) * u.hour
        logger.info(f"Will terminate after {max_wait_hours} hours.")
        self.t_start = Time.now()

        self.email_wait_hours = float(email_wait_hours) * u.hour

        if np.sum(check_email) == 2:

            if self.email_wait_hours > self.max_wait_hours:
                logger.warning(f"Email was set to {self.email_wait_hours}, "
                               f"but the monitor has a shorter termination period of {self.max_wait_hours}. "
                               f"Setting email to 95% of max wait.")
                self.email_wait_hours = 0.95 * self.max_wait_hours

            logger.info(f"Will send an email summary after {self.email_wait_hours} hours.")
            self.email_info = (email_sender, email_recipients)
            self.email_to_send = True

        else:
            logger.info("No email notification configured.")
            self.email_info = None
            self.email_to_send = False

        self.processed_science = []

        # default to "pipeline default cal requirements"

        if cal_requirements is None:
            cal_requirements = self.pipeline.default_cal_requirements

        if cal_requirements is not None:
            self.cal_images, self.cal_headers = find_required_cals(
                latest_dir=str(self.raw_image_directory),
                night=night,
                open_f=self.pipeline.load_raw_image,
                requirements=cal_requirements
            )
        else:
            self.cal_images = self.cal_headers = []

    def summarise_errors(
            self,
            errorstack: ErrorStack,
    ):

        error_summary = errorstack.summarise_error_stack(verbose=False)
        summary = f"Processed a total of {len(self.processed_science)} science images. \n\n" + error_summary + " \n"

        logger.info(f"Writing error log to {self.error_path}")
        errorstack.summarise_error_stack(verbose=True, output_path=self.error_path)

        if self.email_info is not None:

            sender, recipients = self.email_info

            subject = f"{self.pipeline_name}: Summary for night {self.night}"

            send_gmail(
                email_sender=sender,
                email_recipients=recipients,
                email_subject=subject,
                email_text=summary,
                attachments=[self.log_path, self.error_path]
            )
        else:
            print(summary)

    def process_full_night(
            self,
    ):
        batches, errorstack = self.pipeline.reduce_images([[[], []]], catch_all_errors=True)
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

        # Set up a worker thread to process database load
        worker = Thread(target=self.process_load_queue, args=(watchdog_queue,))
        worker.daemon = True
        worker.start()

        # setup watchdog to monitor directory for trigger files
        logger.info(f"Watching {self.raw_image_directory}")

        event_handler = NewImageHandler(watchdog_queue)
        observer = Observer()
        observer.schedule(event_handler, path=str(self.raw_image_directory))
        observer.start()

        try:
            while (Time.now() - self.t_start) < self.max_wait_hours:
                time.sleep(2)
        finally:
            logger.info(f"No longer waiting for new images.")
            observer.stop()
            observer.join()
            self.update_error_log()

    def update_error_log(self):
        self.errorstack.summarise_error_stack(verbose=True, output_path=self.error_path)

    def process_load_queue(self, q):
        '''This is the worker thread function. It is run as a daemon
           threads that only exit when the main thread ends.

           Args
           ==========
             q:  Queue() object
        '''
        while True:

            if self.email_to_send:
                if Time.now() - self.t_start > self.email_wait_hours:
                    logger.info(f"More than {self.email_wait_hours} hours have elapsed. Sending summary email.")
                    self.summarise_errors(errorstack=self.errorstack)
                    self.email_to_send = False

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

                    img, header = self.pipeline.load_raw_image(event.src_path)

                    is_science = header["OBSCLASS"] == "science"

                    all_img = [img] + copy.deepcopy(self.cal_images)
                    all_headers = [header] + copy.deepcopy(self.cal_headers)

                    if not is_science:
                        print(f"Skipping {event.src_path} (calibration image)")
                    else:
                        print(f"Reducing {event.src_path} (science image)")
                        _, errorstack = self.pipeline.reduce_images(
                            batches=[[all_img, all_headers]],
                            selected_configurations=self.realtime_configurations,
                            catch_all_errors=True
                        )
                        self.processed_science.append(event.src_path)
                        self.errorstack += errorstack
                        self.update_error_log()

            else:
                time.sleep(1)
