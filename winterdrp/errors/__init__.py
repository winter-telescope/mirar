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
        self.non_critical_bool = isinstance(self.error, NoncriticalProcessingError)

    def message_known_error(self) -> str:
        return f"This error {['was not', 'was'][self.known_error_bool]} a known error raised by winterdrp."

    def generate_log_message(self) -> str:
        return f"Error for processor {self.processor_name} at time {self.t_error} UT: " \
               f"{type(self.error).__name__} affected batch of length {len(self.contents)}. " \
               f"{self.message_known_error()}. \n"

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
        self.noncritical_reports = []
        self.failed_images = []

        if reports is not None:
            for report in reports:
                self.add_report(report)

    def add_report(self, report: ErrorReport):
        if report.non_critical_bool:
            self.noncritical_reports.append(report)
        else:
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
                  f"{int(len(is_known_error) - np.sum(is_known_error))}/{len(is_known_error)} " \
                  f"errors were errors not raised by winterdrp.\n " \
                  f"The remaining {int(np.sum(is_known_error))}/{len(is_known_error)} errors were known errors " \
                  f"raised by winterdrp. \n" \
                  f"An additional {len(self.noncritical_reports)} non-critical errors were raised. \n" \

        all_reports = self.reports + self.noncritical_reports

        if len(all_reports) > 0:

            summary += f"The following {len(self.failed_images)} images were affected " \
                       f"by at least one error during processing: \n " \
                       f"{self.failed_images} \n \n" \
                       f"Summarising each error: \n\n"

            logger.error(f"Found {len(self.reports)} errors caught by code.")

            errors = [type(x.error).__name__ for x in all_reports]

            for error_type in list(set(errors)):
                summary += f"Found {errors.count(error_type)} counts of error {error_type}. \n"

            if verbose:
                for report in all_reports:
                    summary += str(report.generate_full_traceback())

        else:
            msg = "No raised errors found in processing"
            logger.info(msg)
            summary += f"\n {msg}"

        if output_path is not None:
            logger.error(f"Saving tracebacks of caught errors to {output_path}")
            with open(output_path, "w") as f:
                f.write(summary)

        return summary
