from winterdrp.paths import base_name_key
import traceback
from datetime import datetime
import logging
import numpy as np

logger = logging.getLogger(__name__)


class BaseProcessorError(Exception):
    pass


class ProcessorError(BaseProcessorError):
    pass


class NoncriticalProcessingError(BaseProcessorError):
    pass


class ImageNotFoundError(ProcessorError, FileNotFoundError):
    pass


class ErrorReport:

    def __init__(
            self,
            error,
            processor_name,
            contents: list[str]
    ):
        self.error = error
        self.processor_name = processor_name
        self.contents = contents
        self.t_error = datetime.now()
        self.known_error_bool = isinstance(self.error, BaseProcessorError)

    def message_known_error(self) -> str:
        return f"This error {['was not', 'was'][self.known_error_bool]} a known error raised by winterdrp."

    def generate_log_message(self) -> str:
        return f"Error for processor {self.processor_name} at time {self.t_error} UT: " \
               f"{type(self.error).__name__} affected batch of length {len(self.contents)}. " \
               f"{self.message_known_error()}" \


    def generate_full_traceback(self) -> str:
        msg = f"Error for processor {self.processor_name} at {self.t_error} (local time): \n " \
              f"{''.join(traceback.format_tb(self.error.__traceback__))}" \
              f"{type(self.error).__name__}: {self.error} \n  " \
              f"This error affected the following files: {self.contents} \n" \
              f"{self.message_known_error()} \n \n"
        return msg


class ErrorStack:

    def __init__(
            self,
            reports: list[ErrorReport] = None
    ):
        self.reports = []
        self.failed_images = []

        if reports is not None:
            for report in reports:
                self.add_report(report)

    def add_report(self, report: ErrorReport):
        self.reports.append(report)
        all_failed_images = self.failed_images + report.contents
        self.failed_images = sorted(list(set(all_failed_images)))

    def __add__(self, other):
        for report in other.reports:
            self.add_report(report)
        return self

    def summarise_error_stack(
            self,
            output_path=None,
            verbose: bool = True
    ) -> str:

        is_known_error = [x.known_error_bool for x in self.reports]

        summary = f"Error report summarising {len(self.reports)} errors. \n \n" \
                  f"{np.sum(is_known_error)}/{len(is_known_error)} errors were known errors " \
                  f"raised by winterdrp. \n" \
                  f"The remaining {len(is_known_error) - np.sum(is_known_error)}/{len(is_known_error)} " \
                  f"errors were other errors not raised by winterdrp.\n  \n" \

        if len(self.reports) > 0:

            summary += f"The following images were affected by at least one error during processing: \n " \
                       f"{self.failed_images} \n \n" \
                       f"Summarising each error: \n\n"

            logger.error(f"Found {len(self.reports)} errors caught by code.")

            for report in self.reports:
                if verbose:
                    summary += str(report.generate_full_traceback())
                else:
                    summary += report.generate_log_message()

            if output_path is not None:

                logger.error(f"Saving tracebacks of caught errors to {output_path}")

                with open(output_path, "w") as f:
                    f.write(summary)

        else:
            msg = "No raised errors found in processing"
            logger.info(msg)
            summary += f"\n {msg}"

        return summary
