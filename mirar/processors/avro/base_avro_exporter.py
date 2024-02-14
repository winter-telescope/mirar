"""
Module with classes to make avro alert packets
"""

import logging
import time
from pathlib import Path

import fastavro
from fastavro.schema import load_schema
from fastavro.types import Schema

from mirar.data import SourceBatch, SourceTable
from mirar.paths import get_output_dir
from mirar.processors.base_processor import BaseSourceProcessor

logger = logging.getLogger(__name__)


class BaseAvroExporter(BaseSourceProcessor):
    """Class to generate Avro Packets from a dataframe of candidates.

    Attributes:
        output_sub_dir (str): output data path.
        base_name (str): 4-letter code for telescope.
        save_local (bool): save avro packets to out_sub_dir.
        broadcast (bool): send to brokers at IPAC.
    """

    base_key = "AVRO"

    def __init__(
        self,
        base_name: str,
        avro_schema_path: Path | str,
        output_sub_dir: str = "avro",
        save_local: bool = True,
        broadcast: bool = False,
    ):
        super().__init__()
        self.output_sub_dir = output_sub_dir
        self.base_name = base_name
        self.avro_schema_path = Path(avro_schema_path)
        self.save_local = save_local
        self.broadcast = broadcast

        self.schema = load_schema(self.avro_schema_path)

    def __str__(self) -> str:
        return (
            f"Creates avro packets with '{self.avro_schema_path.name}' schema, "
            f"and save them to '{self.output_sub_dir}' directory."
        )

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        for source_table in batch:
            new_alerts, save_path, topic_name = self.make_alerts(source_table)
            self.process_alerts(new_alerts, save_path, topic_name)
        return batch

    def make_alerts(
        self, source_table: SourceTable
    ) -> tuple[list[dict], Path, str | None]:
        """
        Make avro alerts from a source table

        :param source_table: input source table
        :return: list of avro alerts, path to save avro alerts to, topic name
        """
        raise NotImplementedError

    def get_sub_output_dir(self):
        """Returns path of output subdirectory."""
        return get_output_dir(self.output_sub_dir, self.night_sub_dir)

    @staticmethod
    def save_alert_packets(packets: list[dict], schema: Schema, save_path: Path | str):
        """
        Saves packets to output path.

        :param packets: list of packets to save
        :param schema: schema of the packets
        :param save_path: path to save packets to
        """
        with open(save_path, "wb") as out:
            fastavro.writer(out, schema, packets)

    @staticmethod
    def _send_alert(topicname, records, schema):
        """
        Function to send alert to Kafka broker

        :param topicname: Name of the topic to send to
        :param records: Records to send
        :param schema: Schema of the records
        :return: None
        """

    def broadcast_single_alert_packet(self, packet, schema, topic_name):
        """
        Sends avro-formatted packets to specified topicname using Kafka.
        Modified from https://github.com/dekishalay/pgirdps

        Args:
            packet (dict): candidate data in avro packed dict format.
            schema (dict): schema definition.
            topic_name (str): name of the topic sending to, e.g. ztf_20191221.

        Returns:
            (int): 1 if broadcast successful or -1 if candidate not sent.
        """
        try:
            self._send_alert(topic_name, [packet], schema)
            return 1
        except OSError:
            logger.warning(f"Could not send candid {packet['candidate']['candid']}")
            return -1

    def process_alerts(
        self, alerts: list[dict], save_path: Path, topic_name: str | None = None
    ):
        """
        Top level method to process avro alerts.

        :param alerts: list of avro alerts
        :param save_path: path to save avro alerts to
        :param topic_name: name of the topic sending to, e.g. ztf_20191221

        :return: None
        """
        # Save avro packets to local directory
        save_path.parent.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Saving {len(alerts)} alerts to {save_path}")
        self.save_alert_packets(alerts, self.schema, save_path)

        if self.broadcast:
            t_start = time.time()

            if topic_name is None:
                raise ValueError("topic_name must be specified if broadcast is True")

            logger.debug(f"Broadcasting {len(alerts)} alerts to {topic_name}")

            # Export avro packets to Kafka broker
            successes = 0  # successful kafka producing of a candidate
            for cand in alerts:
                flag = self.broadcast_single_alert_packet(cand, self.schema, topic_name)
                if flag > 0:
                    successes += 1

            t_end = time.time()
            logger.debug(
                f"Took {(t_end - t_start):.2f} seconds to process. "
                f"{successes} of {len(alerts)} successfully broadcast."
            )
