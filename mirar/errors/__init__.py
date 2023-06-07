"""
Central module for handling errors during processing.

In general, the philosophy is that :class:`~mirar.data.base_data.DataBatch` objects
should be processed, unless an error occurs.

If there is an error, then an :class:`~mirar.errors.error_report.ErrorReport`
should be created for the error.

These :class:`~mirar.errors.error_report.ErrorReport` objects should be collated in
a single :class:`~mirar.errors.error_stack.ErrorStack` object.

After processing is complete, the :class:`~mirar.errors.error_stack.ErrorStack` can
then be used to summarise these errors, and track which images failed.

Ideally, it should be understood why the processing failed for a given image. Therefore
the code distinguishes between **internal errors** which were deliberately raised, and
**external errors** which were not deliberately raised.

.. include:: ../../mirar/errors/exceptions.py
    :start-line: 3
    :end-line: 14

"""
from mirar.errors.error_report import ErrorReport
from mirar.errors.error_stack import ErrorStack
from mirar.errors.exceptions import (
    ImageNotFoundError,
    NoncriticalProcessingError,
    ProcessorError,
)
