import logging

import numpy as np

from winterdrp.errors.error_report import ErrorReport
from winterdrp.paths import __version__, package_name

logger = logging.getLogger(__name__)


class ErrorStack:
    def __init__(self, reports: list[ErrorReport] = None):
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

    def get_all_reports(self) -> list[ErrorReport]:
        return self.reports + self.noncritical_reports

    def summarise_error_stack(self, output_path=None, verbose: bool = True) -> str:

        is_known_error = [x.known_error_bool for x in self.reports]

        summary = (
            f"Error report summarising {len(self.reports)} errors. \n"
            f"Code version: {package_name}=={__version__} \n \n"
            f"{int(len(is_known_error) - np.sum(is_known_error))}/{len(is_known_error)} "
            f"errors were errors not raised by {package_name}. \n"
            f"The remaining {int(np.sum(is_known_error))}/{len(is_known_error)} errors were known errors "
            f"raised by {package_name}. \n"
            f"An additional {len(self.noncritical_reports)} non-critical errors were raised. \n"
        )
        all_reports = self.get_all_reports()

        if len(all_reports) > 0:

            summary += (
                f"The following {len(self.failed_images)} images were affected "
                f"by at least one error during processing: \n "
            )

            if verbose:
                summary += f"{self.failed_images} \n \n"

            summary += f"Summarising each error: \n\n"

            logger.error(f"Found {len(self.reports)} errors caught by code.")

            error_lines = [x.get_error_message() for x in all_reports]

            for error_type in list(set(error_lines)):

                matching_errors = [
                    x for x in all_reports if x.get_error_message() == error_type
                ]

                img_paths = []
                error_name = None
                for x in matching_errors:
                    img_paths += x.contents
                    if error_name is None:
                        error_name = x.get_error_name()
                img_paths = list(set(img_paths))

                line = error_type.split("\n")[0]
                summary += (
                    f"Found {error_lines.count(error_type)} counts of error {error_name}, "
                    f"affecting {len(img_paths)} images: \n{line}.\n \n"
                )

            summary += " \n"

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
