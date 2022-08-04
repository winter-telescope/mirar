import subprocess
import os
import datetime
from watchdog.events import FileSystemEventHandler
from threading import Thread
from queue import Queue
import time
import sys
from watchdog.observers import Observer
from winterdrp.pipelines import get_pipeline
from winterdrp.errors import ErrorStack
from winterdrp.utils.send_email import send_gmail
from winterdrp.paths import get_output_path, raw_img_dir
import numpy as np
import logging
from astropy.time import Time
from astropy import units as u
from winterdrp.pipelines.summer.summer_pipeline import load_raw_summer_image
from winterdrp.processors.utils.image_selector import select_from_images
from winterdrp.processors.utils.image_loader import load_from_dir, ImageNotFoundError
from winterdrp.processors.utils.supplement_cals import CalHunter, CalRequirement

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
            log_level="INFO",
    ):

        self.errorstack = ErrorStack
        self.night = night
        self.pipeline_name = pipeline

        if not isinstance(realtime_configurations, list):
            realtime_configurations = [realtime_configurations]
        self.realtime_configurations = realtime_configurations

        self.pipeline = get_pipeline(pipeline, night=night, selected_configurations=realtime_configurations)

        self.raw_image_directory = raw_img_dir(sub_dir=self.pipeline.night_sub_dir)

        self.log_level = log_level
        self.log_path = self.configure_logs(log_level)

        self.processed_science = []

        if cal_requirements is not None:
            self.cal_images, self.cal_headers = self.get_latest_cals(cal_requirements)
        else:
            self.cal_images = self.cal_headers = []

        # self.cal_cache = self.get_latest_cals(required_cals)

        check_email = np.sum([x is not None for x in [email_recipients, email_sender]])
        if np.sum(check_email) == 1:
            err = "In order to send emails, you must specify both a a sender and a recipient. \n" \
                  f"In this case, sender is {email_sender} and recipent is {email_recipients}."
            logger.error(err)
            raise ValueError(err)

        if np.sum(check_email) == 2:
            self.email_info = (email_sender, email_recipients)
        else:
            self.email_info = None

    def summarise_errors(
            self,
            errorstack: ErrorStack,
    ):

        summary = errorstack.summarise_error_stack()

        if self.email_info is not None:

            sender, recipients = self.email_info

            subject = f"{self.pipeline_name}: Summary for night {self.night}"

            send_gmail(
                email_sender=sender,
                email_recipients=recipients,
                email_subject=subject,
                email_text=summary,
                # attachments=[self.log_path]
            )
        else:
            print(summary)

    def process_full_night(
            self,
    ):
        batches, errorstack = self.pipeline.reduce_images([[[], []]], catch_all_errors=True)
        self.summarise_errors(errorstack=errorstack)

    def configure_logs(self, log_level="INFO"):

        log_output_path = get_output_path(
            base_name=f"{self.night}_processing_log.txt",
            dir_root=self.pipeline.night_sub_dir,
        )

        try:
            os.makedirs(os.path.dirname(log_output_path))
        except OSError:
            pass

        # log = logging.getLogger("winterdrp")
        #
        # handler = logging.FileHandler(log_output_path)
        # # handler = logging.StreamHandler(sys.stdout)
        # formatter = logging.Formatter('%(asctime)s: %(name)s [l %(lineno)d] - %(levelname)s - %(message)s')
        # handler.setFormatter(formatter)
        # log.addHandler(handler)
        # log.setLevel(log_level)

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
        # event_handler = SqlLoaderWatchdog(watchdog_queue)

        print(f"Watching {self.raw_image_directory}")

        event_handler = NewImageHandler(watchdog_queue)
        observer = Observer()
        observer.schedule(event_handler, self.raw_image_directory, recursive=True)
        observer.start()

        try:
            while True:
                time.sleep(2)
        except KeyboardInterrupt:
            observer.stop()

        observer.join()

    def process_load_queue(self, q):
        '''This is the worker thread function. It is run as a daemon
           threads that only exit when the main thread ends.

           Args
           ==========
             q:  Queue() object
        '''
        while True:
            if not q.empty():
                event = q.get()
                now = datetime.datetime.utcnow()

                print(now)

                img, header = load_raw_summer_image(event.src_path)

                is_science = header["OBSCLASS"] == "science"

                all_img = [img] + self.cal_images
                all_headers = [header] + self.cal_headers

                if not is_science:
                    pass
                    # self.raw_cals.append(event.src_path)
                else:
                    _, errorstack = self.pipeline.reduce_images(
                        batches=[[all_img, all_headers]],
                        selected_configurations=self.realtime_configurations,
                        catch_all_errors=False
                    )
                    self.processed_science.append(event.src_path)
                    self.errorstack += errorstack
            else:
                time.sleep(1)

    def get_latest_cals(self, requirements):

        latest_dir = self.raw_image_directory

        preceding_dirs = []

        for x in os.listdir(os.path.dirname(os.path.dirname(latest_dir))):
            if x[0] not in ["."]:
                if len(str(x)) == len(str(self.night)):
                    try:
                        if not float(x) > float(self.night):
                            preceding_dirs.append(x)
                    except ValueError:
                        pass

        ordered_nights = sorted(preceding_dirs)[::-1]

        while np.sum([x.success for x in requirements]) != len(requirements):

            if len(ordered_nights) == 0:
                raise ImageNotFoundError("Ran out of nights!")

            new_latest_night = ordered_nights[0]
            ordered_nights = ordered_nights[1:]

            try:

                logger.info(f"Checking night {new_latest_night}")

                dir_to_load = self.raw_image_directory.replace(self.night, new_latest_night)

                images, headers = load_from_dir(
                    dir_to_load, open_f=load_raw_summer_image
                )

                for requirement in requirements:
                    if not requirement.success:
                        requirement.check_images(images, headers)

            except ImageNotFoundError:
                pass

        all_images = []
        all_headers = []

        for requirement in requirements:
            for key, (imgs, headers) in requirement.data.items():
                print(key, len(imgs))
                all_images += imgs
                all_headers += headers

        logger.info(f"Found {len(all_images)} calibration images")

        return all_images, all_headers


if __name__ == '__main__':


    # root = logging.getLogger()
    # root.setLevel(logging.DEBUG)
    #
    # handler = logging.StreamHandler(sys.stdout)
    # handler.setLevel(logging.DEBUG)
    # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # handler.setFormatter(formatter)
    # root.addHandler(handler)

    ln = Time.now() - 1. * u.day
    last_night = str(ln).split(" ")[0].replace("-", "")

    required_cals = [
        CalRequirement(target_name="bias", required_field="EXPTIME", required_values=["0.0"]),
        CalRequirement(target_name="flat", required_field="FILTERID", required_values=["u", "g", "r", "i"]),
    ]

    scrutineer = Monitor(
        pipeline="summer",
        night=last_night,
        cal_requirements=required_cals,
        realtime_configurations=["realtime"],
    )
    scrutineer.process_realtime()
