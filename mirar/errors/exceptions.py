"""
Module containing common exceptions or base exceptions for the code.

In general, all internal errors should inherit from the
:class:`mirar.errors.exceptions.BaseProcessorError` class.

If the error is critical (i.e the image should not be processed further), then an
error should be raised which inherits from the
:class:`mirar.errors.exceptions.ProcessorError` class.

If the error is non-critical (so processing should continue), then an
error should be raised which inherits from the
:class:`mirar.errors.exceptions.NoncriticalProcessingError` class. In that case,
processing will continue but the error will be logged.
"""


class BaseProcessorError(Exception):
    """
    Base exception, from which all internal exceptions should derive
    """


class ProcessorError(BaseProcessorError):
    """
    Base class for all critical internal exceptions
    """


class NoncriticalProcessingError(BaseProcessorError):
    """
    Base class for all non-critical internal exceptions
    """


class ImageNotFoundError(ProcessorError, FileNotFoundError):
    """
    Base class for all exceptions concerning missing images
    """
