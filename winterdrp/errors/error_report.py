import logging
import traceback
from datetime import datetime

from winterdrp.errors.exceptions import BaseProcessorError, NoncriticalProcessingError

logger = logging.getLogger(__name__)


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
               f"{self.get_error_name()} affected batch of length {len(self.contents)}. " \
               f"{self.message_known_error()}. \n"

    def generate_full_traceback(self) -> str:
        msg = f"Error for processor {self.processor_name} at {self.t_error} (local time): \n " \
              f"{''.join(traceback.format_tb(self.error.__traceback__))}" \
              f"{self.get_error_name()}: {self.error} \n  " \
              f"This error affected the following files: {self.contents} \n" \
              f"{self.message_known_error()} \n \n"
        return msg

    def get_error_name(self) -> str:
        return type(self.error).__name__

    def get_error_message(self) -> str:
        return traceback.format_tb(self.error.__traceback__)[-1]

    def get_error_line(self) -> str:
        return self.get_error_message().split('\n')[0]