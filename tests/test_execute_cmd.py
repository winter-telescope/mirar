"""
Module for testing the execution of external commands
"""

import logging

from mirar.testing import BaseTestCase
from mirar.utils.execute_cmd import ExecutionError, decode_output, run_local

logger = logging.getLogger(__name__)

ERROR_MESSAGE = "this went wrong"


class TestExecuteCmd(BaseTestCase):
    """Class for testing ..func::run_local"""

    def setUp(self):
        self.logger = logging.getLogger(__name__)

    def test_captures_output(self):
        """
        Test that a successful command writes nothing to the terminal

        The external tools are far more verbose than mirar itself, so both of
        their streams are captured rather than inherited.

        :return: None
        """
        self.logger.info("\n\n Testing command execution \n\n")

        with self.assertLogs("mirar.utils.execute_cmd", level="DEBUG") as logs:
            run_local(f"echo '{ERROR_MESSAGE}' >&2; echo fine")

        output = "".join(logs.output)

        self.assertIn("fine", output)
        self.assertIn(ERROR_MESSAGE, output)

    def test_reports_error_output(self):
        """
        Test that the output of a failed command reaches the error

        :return: None
        """
        with self.assertRaises(ExecutionError) as err:
            run_local(f"echo '{ERROR_MESSAGE}' >&2; exit 1")

        self.assertIn(ERROR_MESSAGE, str(err.exception))

    def test_decode_output(self):
        """
        Test that output is decoded whether or not it was captured

        :return: None
        """
        self.assertEqual(decode_output(b"bytes"), "bytes")
        self.assertEqual(decode_output("text"), "text")
        self.assertIn("No output", decode_output(None))
