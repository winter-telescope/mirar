"""
Utility for retrying flaky operations (e.g. external network calls) with
exponential backoff.
"""

import functools
import logging
import time
from typing import Callable, TypeVar

logger = logging.getLogger(__name__)

T = TypeVar("T")  # pylint: disable=invalid-name

DEFAULT_RETRY_ATTEMPTS = 5
DEFAULT_RETRY_BACKOFF_BASE = 2.0


def retry_on_exception(
    exceptions: tuple[type[Exception], ...] = (Exception,),
    attempts: int = DEFAULT_RETRY_ATTEMPTS,
    backoff_base: float = DEFAULT_RETRY_BACKOFF_BASE,
) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """
    Decorator to retry a function with exponential backoff, if it raises
    one of `exceptions`. Useful for flaky external queries (e.g. an
    astronomical catalog service), which can fail transiently.

    :param exceptions: Exception types that should trigger a retry
    :param attempts: Maximum number of attempts
    :param backoff_base: Base for the exponential backoff delay
        (backoff_base ** attempt seconds)
    :return: Decorated function
    """

    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @functools.wraps(func)
        def wrapper(*args, **kwargs) -> T:
            for attempt in range(attempts):
                try:
                    return func(*args, **kwargs)
                except exceptions as exc:
                    if attempt == attempts - 1:
                        logger.error(
                            f"{func.__qualname__} failed after {attempts} attempts"
                        )
                        raise
                    delay = backoff_base**attempt
                    logger.warning(
                        f"{func.__qualname__} failed with {exc!r} "
                        f"(attempt {attempt + 1}/{attempts}), "
                        f"retrying in {delay:.1f}s"
                    )
                    time.sleep(delay)
            raise AssertionError("unreachable")  # pragma: no cover

        return wrapper

    return decorator
