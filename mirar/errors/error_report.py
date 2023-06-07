"""
Module for ErrorReport objects.

An :class:`~mirar.errors.error_report.ErrorReport` object summarises a single error
raised by the code.
"""
import logging
import traceback
from datetime import datetime
from pathlib import Path

from mirar.errors.exceptions import BaseProcessorError, NoncriticalProcessingError

logger = logging.getLogger(__name__)


class ErrorReport:
    """
    Class representing a single error raised during processing
    """

    def __init__(
        self, error: Exception, processor_name: str, contents: list[str] | list[Path]
    ):
        self.error = error
        self.processor_name = processor_name
        self.contents = contents
        self.t_error = datetime.now()
        self.known_error_bool = isinstance(self.error, BaseProcessorError)
        self.non_critical_bool = isinstance(self.error, NoncriticalProcessingError)

    def message_known_error(self) -> str:
        """
        Returns a human-readable string describing whether the error was an internal one
        intentionally raised by the code, or an unexpected external error

        :return: String describing whether error was known
        """
        return (
            f"This error {['was not', 'was'][self.known_error_bool]} "
            f"a known error raised by mirar."
        )

    def generate_log_message(self) -> str:
        """
        Returns a human-readable string describing high-livel details about the error.

        :return: String summary
        """
        return (
            f"Error for processor {self.processor_name} at time {self.t_error} UT: "
            f"{self.get_error_name()} affected batch of length {len(self.contents)}. "
            f"{self.message_known_error()}. \n"
        )

    def generate_full_traceback(self) -> str:
        """
        Returns a verbose string summarising the error

        :return: String
        """
        msg = (
            f"Error for processor {self.processor_name} at {self.t_error} "
            f"(local time): \n "
            f"{''.join(traceback.format_tb(self.error.__traceback__))}"
            f"{self.get_error_name()}: {self.error} \n  "
            f"This error affected the following files: {self.contents} \n"
            f"{self.message_known_error()} \n \n"
        )
        return msg

    def get_error_name(self) -> str:
        """
        Returns the name of the error

        :return: Name
        """
        return type(self.error).__name__

    def get_error_message(self) -> str:
        """
        Returns the full error message raised by python

        :return: String for single line
        """
        return traceback.format_tb(self.error.__traceback__)[-1]

    def get_error_line(self) -> str:
        """
        Returns only the critical error line raised by python

        :return: string
        """
        return self.get_error_message().split("\n")[0]
