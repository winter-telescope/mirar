"""
Module with classes to make avro alert packets
"""

import io
import logging
from pathlib import Path

import confluent_kafka
import fastavro
import pandas as pd
from fastavro.types import Schema

from mirar.data import SourceTable
from mirar.paths import BASE_NAME_KEY
from mirar.processors.avro.base_avro_exporter import BaseAvroExporter

logger = logging.getLogger(__name__)


class IPACAvroExporter(BaseAvroExporter):
    """Class to generate Avro Packets from a dataframe of candidates.

    Attributes:
        output_sub_dir (str): output data path.
        base_name (str): 4-letter code for telescope.
        save_local (bool): save avro packets to out_sub_dir.
        broadcast (bool): send to brokers at IPAC.
    """

    def __init__(self, *args, topic_prefix: str | None = None, **kwargs):
        super().__init__(*args, **kwargs)
        self.alert_schema, self.candidate_schema, self.prv_schema = self._load_schemas()
        self.topic_prefix = topic_prefix

    def _load_schemas(self) -> tuple[Schema, Schema, Schema]:
        """
        Unpack an IPAC-style avro schema (alert, alert.candidate,
        & alert.prv_candidate) its three components

        :return: alert, candidate, prv_candidate schemas
        """
        keys = sorted(self.schema["__named_schemas"].keys())

        assert len(keys) == 3, "Unrecognised alert structure"

        alert_schema = self.schema["__named_schemas"][keys[0]]
        assert alert_schema["name"].split(".")[-1] == "alert", ""
        candidate_schema = self.schema["__named_schemas"][keys[1]]
        assert candidate_schema["name"].split(".")[-1] == "candidate"
        prv_candidate_schema = self.schema["__named_schemas"][keys[2]]
        assert prv_candidate_schema["name"].split(".")[-1] == "prv_candidate"

        return alert_schema, candidate_schema, prv_candidate_schema

    @staticmethod
    def fill_schema(schema: Schema, row: pd.Series, metadata: dict) -> dict:
        """
        Fill an avro schema with data from a row of a pandas dataframe

        :param schema: Schema to fill
        :param row: Row of pandas dataframe
        :param metadata: Metadata to fill
        :return: Dictionary of filled schema
        """
        new = {}

        for field in schema["fields"]:
            key = field["name"]
            if key in row.keys():
                new[key] = row[key]
            elif key.upper() in row.keys():
                new[key] = row[key.upper()]
            elif key in metadata.keys():
                new[key] = metadata[key]
            elif key.upper() in metadata.keys():
                new[key] = metadata[key.upper()]

        return new

    def make_alerts(
        self, source_table: SourceTable
    ) -> tuple[list[dict], Path, str | None]:
        """
        Make avro alerts from a source table

        :param source_table: input source table
        :return: list of avro alerts
        """
        new_alerts = []

        metadata = source_table.get_metadata()

        for _, row in source_table.get_data().iterrows():
            alert = self.fill_schema(self.alert_schema, row, metadata)
            candidate = self.fill_schema(self.candidate_schema, row, metadata)
            alert["candidate"] = candidate

            prv_candidates = []
            if "prv_candidates" in row.keys():
                prv_cands = pd.DataFrame(row["prv_candidates"])
                if len(prv_candidates) > 0:
                    for _, prv_row in prv_cands.iterrows():
                        prv_dict = {}
                        for key in self.prv_schema["fields"]:
                            if key in prv_row.keys():
                                prv_dict[key] = prv_row[key]
                        prv_candidates.append(prv_dict)
            alert["prv_candidate"] = prv_candidates

            new_alerts.append(alert)

        save_path = self.get_sub_output_dir().joinpath(
            Path(metadata[BASE_NAME_KEY]).with_suffix(".avro")
        )

        topic_name = self.get_topic_name()

        return new_alerts, save_path, topic_name

    def get_topic_name(self) -> str:
        """
        Get the topic name for a source table

        :return: topic name
        """
        if self.topic_prefix is not None:
            topic_name = f"{self.topic_prefix}_{self.night}"
        else:
            topic_name = None
        return topic_name

    @staticmethod
    def _send_alert(topicname, records, schema):
        """Send an avro "packet" to a particular topic at IPAC.
        Modified from: https://github.com/dekishalay/pgirdps

        Args:
            topicname (str): name of the topic sending to,
            e.g. ztf_20191221_programid2_zuds.
            records (list): a list of dictionaries (the avro packet to send).
            schema (dict): schema definition.
        """
        with io.BytesIO() as out:
            fastavro.writer(out, schema, records)
            out.seek(0)  # go back to the beginning

            # Connect to the IPAC Kafka brokers
            producer = confluent_kafka.Producer(
                {
                    "bootstrap.servers": "ztfalerts04.ipac.caltech.edu:9092,"
                    "ztfalerts05.ipac.caltech.edu:9092,"
                    "ztfalerts06.ipac.caltech.edu:9092"
                }
            )

            # Send an avro alert
            producer.produce(topic=topicname, value=out.read())
            producer.flush()
