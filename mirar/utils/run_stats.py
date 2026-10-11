"""
Module for RunStats objects.

A :class:`~mirar.utils.run_stats.RunStats` object records how long each step of a
pipeline run took, together with the circumstances of the run, so that runs can be
compared against each other.

The statistics are snapshotted as each step finishes, rather than read back from the
processors afterwards. Processor objects are shared between runs of the same pipeline,
so their attributes only ever describe the most recent run of that processor.
"""

import logging
import socket
import subprocess
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import pandas as pd

from mirar.paths import __version__

if TYPE_CHECKING:
    from mirar.processors.base_processor import BaseProcessor

logger = logging.getLogger(__name__)

GIT_TIMEOUT = 10

# Timings are recorded to the millisecond, which is well beyond the run-to-run noise
RUNTIME_PRECISION = 3

RUNTIME_COLUMNS = [
    "runtime",
    "batch_runtime_total",
    "batch_runtime_max",
    "runtime_per_block",
]


def get_git_commit() -> str:
    """
    Get the git commit of the code being run, if it is a git checkout

    An uncommitted working tree is flagged as '-dirty', because the package version
    cannot distinguish two working trees built from the same release. Untracked files
    are ignored, since they are not part of the code being run.

    :return: Short commit hash, or an empty string if this is not a git checkout
    """
    repo_dir = Path(__file__).parents[2]

    git_args = ["git", "-C", str(repo_dir)]

    try:
        commit = subprocess.run(
            git_args + ["rev-parse", "--short", "HEAD"],
            capture_output=True,
            check=True,
            text=True,
            timeout=GIT_TIMEOUT,
        ).stdout.strip()

        modified = subprocess.run(
            git_args + ["status", "--porcelain", "--untracked-files=no"],
            capture_output=True,
            check=True,
            text=True,
            timeout=GIT_TIMEOUT,
        ).stdout.strip()

    except (subprocess.SubprocessError, OSError):
        logger.debug(f"Could not get a git commit for {repo_dir}")
        return ""

    return f"{commit}-dirty" if modified else commit


class RunStats:
    """
    Container class to hold the statistics of each step of a single pipeline run
    """

    def __init__(self, instrument: str, night: str):
        self.instrument = instrument
        self.night = night
        self.version = __version__
        self.commit = get_git_commit()
        self.hostname = socket.gethostname()
        self.t_start = datetime.now()
        self.rows: list[dict] = []

    def add_step(self, configuration: str, step: int, processor: "BaseProcessor"):
        """
        Record the statistics of a step which has just finished

        :param configuration: Name of the pipeline configuration being run
        :param step: Position of the step within its configuration, counting from one
        :param processor: Processor which has just been applied
        :return: None
        """
        self.rows.append(
            {
                "instrument": self.instrument,
                "night": self.night,
                "version": self.version,
                "commit": self.commit,
                "hostname": self.hostname,
                "started": self.t_start.isoformat(timespec="seconds"),
                "configuration": configuration,
                "step": step,
            }
            | processor.latest_stats()
        )

    @property
    def total_runtime(self) -> float:
        """
        Total time spent in the steps which have been recorded so far

        :return: Runtime in seconds
        """
        return sum(row["runtime"] for row in self.rows)

    def to_dataframe(self) -> pd.DataFrame:
        """
        Build a table of the recorded steps, one row per step

        :return: Dataframe of run statistics
        """
        stats_df = pd.DataFrame(self.rows)

        if len(stats_df) > 0:
            stats_df[RUNTIME_COLUMNS] = stats_df[RUNTIME_COLUMNS].round(
                RUNTIME_PRECISION
            )

        return stats_df

    def summarise_run_stats(self, output_path: Optional[Path] = None) -> pd.DataFrame:
        """
        Summarise the run, and optionally write the statistics as a tsv

        :param output_path: Path to write the statistics to, in .tsv format (optional)
        :return: Dataframe of run statistics
        """
        stats_df = self.to_dataframe()

        if len(self.rows) > 0:
            slowest = max(self.rows, key=lambda row: row["runtime"])
            logger.info(
                f"Ran {len(stats_df)} steps in {self.total_runtime:.1f}s, "
                f"of which the slowest was {slowest['processor']} "
                f"at {slowest['runtime']:.1f}s"
            )

        if output_path is not None:
            logger.info(f"Saving run statistics to {output_path}")
            stats_df.to_csv(output_path, sep="\t", index=False)

        return stats_df
