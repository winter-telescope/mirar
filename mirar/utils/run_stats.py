"""
Module for RunStats objects.

A :class:`~mirar.utils.run_stats.RunStats` object records how long each step of a
pipeline run took, together with the circumstances of the run, so that runs can be
compared against each other.

The results are written as a pair of files sharing a stem: a tsv holding one row per
step, and a json holding the values which describe the run as a whole. Keeping the
run-level values out of the tsv stops them being repeated on every row.
:func:`~mirar.utils.run_stats.read_run_stats` reads the pair back as a single table.

The statistics are snapshotted as each step finishes, rather than read back from the
processors afterwards. Processor objects are shared between runs of the same pipeline,
so their attributes only ever describe the most recent run of that processor.
"""

import json
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

# Rounded as they are recorded, so that a total always matches the rows it sums
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
        stats = processor.latest_stats()

        for key in RUNTIME_COLUMNS:
            stats[key] = round(stats[key], RUNTIME_PRECISION)

        self.rows.append({"configuration": configuration, "step": step} | stats)

    @property
    def metadata(self) -> dict:
        """
        The values describing the run as a whole, which every step shares

        :return: Dictionary of run-level values
        """
        return {
            "instrument": self.instrument,
            "night": self.night,
            "version": self.version,
            "commit": self.commit,
            "hostname": self.hostname,
            "started": self.t_start.isoformat(timespec="seconds"),
            "total_runtime": round(self.total_runtime, RUNTIME_PRECISION),
        }

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
        return pd.DataFrame(self.rows)

    def summarise_run_stats(self, output_path: Optional[Path] = None) -> pd.DataFrame:
        """
        Summarise the run, and optionally write the statistics out

        Writes the per-step table to the given path, and the run-level values to a
        json file beside it, sharing the same stem.

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

            run_path = output_path.with_suffix(".json")

            with open(run_path, "w", encoding="utf-8") as run_file:
                json.dump(self.metadata, run_file, indent=4)

        return stats_df


def read_run_stats(stats_paths: Path | list[Path]) -> pd.DataFrame:
    """
    Read the statistics of one or more runs back as a single table

    The run-level values are stored once per run, and are added back here as columns,
    so that steps from different runs can be told apart after they are combined.

    :param stats_paths: Path, or paths, of the tsv files to read
    :return: Dataframe of run statistics
    """
    if isinstance(stats_paths, Path):
        stats_paths = [stats_paths]

    all_stats = []

    for stats_path in stats_paths:
        stats_df = pd.read_csv(stats_path, sep="\t")

        with open(stats_path.with_suffix(".json"), "r", encoding="utf-8") as run_file:
            metadata = json.load(run_file)

        # Inserted in reverse, so that they lead the table in their original order
        for key, value in reversed(list(metadata.items())):
            stats_df.insert(0, key, value)

        all_stats.append(stats_df)

    return pd.concat(all_stats, ignore_index=True)
