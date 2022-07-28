from winterdrp.paths import base_name_key
import traceback
from astropy.time import Time


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
        msg = f"Error for processor {self.processor_name} at time {self.t_error} UT: \n " \
              f"{traceback.print_exception(self.error)} \n " \
              f"This error affected the following files: {self.contents} \n"
        return msg


def summarise_error_reports(
        reports: list[ErrorReport],
        output_path=None
) -> str:
    summary = "",

    for report in reports:
        summary += report.generate_full_traceback()

    if output_path is not None:
        with open(output_path, "wb") as f:
            f.write(summary)

    return summary
