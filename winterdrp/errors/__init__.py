from winterdrp.paths import base_name_key
import traceback
from astropy.time import Time
import logging

logger = logging.getLogger(__name__)


class ProcessorError(BaseException):
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
        self.t_error = Time.now()

    def generate_log_message(self):
        return f"Error for processor {self.processor_name} at time {self.t_error} UT: " \
               f"{type(self.error).__name__} affected batch of length {len(self.contents)}."

    def generate_full_traceback(self):

        # print(vars(self.error))
        # print(self.error)
        # print(dir(self.error))
        # # print(self.error.with_traceback())
        # print(self.error.__traceback__)
        # print(traceback.format_tb(self.error.__traceback__))
        # print("um")

        msg = f"Error for processor {self.processor_name} at time {self.t_error} UT: \n " \
              f"{''.join(traceback.format_tb(self.error.__traceback__))} \n " \
              f"This error affected the following files: {self.contents} \n"
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
            output_path=None
    ) -> str:

        summary = f"Error report summarising {len(self.reports)} errors. \n" \

        if len(self.reports) > 0:

            summary += f"The following images were affected by at least one error during processing \n: " \
                       f"{self.failed_images} \n \n" \
                       f"Summarising each error: \n"

            logger.error(f"Found {len(self.reports)} errors caught by code.")

            for report in self.reports:
                summary += str(report.generate_full_traceback())

            if output_path is not None:

                logger.error(f"Saving tracebacks of caught errors to {output_path}")

                with open(output_path, "w") as f:
                    f.write(summary)

        else:
            msg = "No raised errors found in processing"
            logger.info(msg)
            summary += f"\n {msg}"

        return summary
