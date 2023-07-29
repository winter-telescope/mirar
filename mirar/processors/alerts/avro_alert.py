"""
Module with classes to make avro alert packets
"""
import datetime
import io
import logging
import time
from pathlib import Path

import avro
import avro.io
import avro.schema
import avro.tool
import confluent_kafka
import fastavro
import pandas as pd
from fastavro.schema import load_schema
from fastavro.types import Schema

from mirar.data import SourceBatch
from mirar.paths import get_output_dir
from mirar.processors.base_processor import BaseSourceProcessor

logger = logging.getLogger(__name__)


class AvroPacketMaker(BaseSourceProcessor):
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
        for source in batch:
            candidate_table = source.get_data()
            self.make_alert(candidate_table)
        return batch

    @staticmethod
    def read_input_df(candidate_df: pd.DataFrame):
        """Takes a DataFrame, which has multiple candidate
        and creates list of dictionaries, each dictionary
        representing a single candidate.

        Args:
            candidate_df (pandas.core.frame.DataFrame): dataframe of all candidates.

        Returns:
            (list[dict]): list of dictionaries, each a candidate.
        """
        all_candidates = []

        for i in range(0, len(candidate_df)):
            candidate = {}
            for key in candidate_df.keys():
                try:
                    if isinstance(candidate_df.iloc[i].get(key), (str, list)):
                        candidate[key] = candidate_df.iloc[i].get(key)
                    else:
                        # change to native python type
                        candidate[key] = candidate_df.iloc[i].get(key).item()
                except AttributeError:  # for IOBytes objs
                    candidate[key] = candidate_df.iloc[i].get(key)
            all_candidates.append(candidate)

        return all_candidates

    @staticmethod
    def write_avro_data(json_dict, avro_schema):
        """Encode json into avro format given a schema.
        For testing packet integrity.

        Args:
            json_dict (dict): data to put into schema, can be dict.
            avro_schema (avro.schema.RecordSchema): output schema.

        Returns:
            (io.BytesIO object): input json in schema format.
        """
        writer = avro.io.DatumWriter(avro_schema)
        bytes_io = io.BytesIO()
        encoder = avro.io.BinaryEncoder(bytes_io)
        writer.write(json_dict, encoder)
        return bytes_io

    @staticmethod
    def read_avro_data(bytes_io, avro_schema):
        """Read avro data and decode with a given schema.
        For testing packet integrity.

        Args:
            bytes_io (_io.BytesIO object): data to put into schema, can be dict.
            avro_schema (avro.schema.RecordSchema): schema to decode from.

        Returns:
            (dict): schema in dictionary format.
        """
        raw_bytes = bytes_io.getvalue()
        bytes_reader = io.BytesIO(raw_bytes)
        decoder = avro.io.BinaryDecoder(bytes_reader)
        reader = avro.io.DatumReader(avro_schema)
        message = reader.read(decoder)
        return message

    @staticmethod
    def create_alert_packet(cand, scicut, refcut, diffcut):
        """Create top level avro packet from input candidate.
        Args:
            cand (dict): all data of a single candidate.
            scicut (bytes): science image cutout.
            refcut (bytes): reference image cutout.
            diffcut (bytes): difference image cutout.

        Returns:
            (dict): schema in dictionary format.
        """
        prev_cands = cand["prv_candidates"]
        alert = {
            "schemavsn": "0.1",
            "publisher": "winter_test",
            "cutoutScience": scicut,
            "cutoutTemplate": refcut,
            "cutoutDifference": diffcut,
            "objectId": cand["objectId"],
            "candid": cand["candid"],
            "candidate": cand,
            "prv_candidates": prev_cands,
        }
        return alert

    def get_sub_output_dir(self):
        """Returns path of output subdirectory."""
        return get_output_dir(self.output_sub_dir, self.night_sub_dir)

    def _save_local(self, candid, records, schema: Schema):
        """
        Save avro file in given schema to output subdirectory.

        Args:
            candid (number type): candidate id.
            records (list): a list of dictionaries
            schema: schema definition.
        """
        avro_packet_path = self.get_sub_output_dir().joinpath(str(candid) + ".avro")
        with open(avro_packet_path, "wb") as out:
            fastavro.writer(out, schema, records)

    def save_alert_packet(self, packet, cand, schema: Schema):
        """Saves packet as .avro to output subdirectory.

        Args:
            packet (dict): candidate data in avro packed dict format.
            cand (dict): all data of a single candidate.
            schema (avro.schema.RecordSchema): schema definition.

        Returns:
            (int): 1 if save successful. Else, -1.
        """
        self._save_local(cand["candid"], [packet], schema)

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
        out = io.BytesIO()

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

    def broadcast_alert_packet(self, packet, cand, schema, topic_name):
        """
        Sends avro-formatted packets to specified topicname using Kafka.
        Modified from https://github.com/dekishalay/pgirdps

        Args:
            packet (dict): candidate data in avro packed dict format.
            cand (dict): all data of a single candidate.
            schema (dict): schema definition.
            topic_name (str): name of the topic sending to,
            e.g. ztf_20191221_programid2_zuds.
            cand_num (int): number of current candidate being sent.
            num_cands (int): total number of candidates to send.

        Returns:
            (int): 1 if broadcast successful or -1 if candidate not sent.
        """
        try:
            self._send_alert(topic_name, [packet], schema)
            return 1
        except OSError:
            logger.info(f"Could not send candid {cand['candid']}")
            return -1

    def _make_ind_packet(self, cand, schema: Schema, topic_name: str):
        """Makes a single avro alert. Saves or broadcasts alert based on global var.

        Args:
            cand (dict): all data of a single candidate.
            schema (dict): schema definition.
            topic_name (str): name of the topic sending to,
            e.g. ztf_20191221_programid2_zuds.

        Returns:
            (int): 1 if broadcast successful,
                   0 if alert neither saved nor broadcast,
                  -1 if candidate not sent or saved.
        """
        # Cutouts are include in the top level alert schema
        scicut = cand.pop("cutoutScience")
        refcut = cand.pop("cutoutTemplate")
        diffcut = cand.pop("cutoutDifference")

        packet = self.create_alert_packet(cand, scicut, refcut, diffcut)

        self.save_alert_packet(packet, cand, schema)

        flag = 1
        if self.broadcast:
            flag += self.broadcast_alert_packet(packet, cand, schema, topic_name)

        return flag

    def _success_message(self):
        """Text for successful save and/or broadcast."""
        message = ""
        if self.save_local:
            message += "saved"
            if self.broadcast:
                message += " and broadcasted"
            return message

        if self.broadcast:
            message += "broadcasted"
        return message

    def make_alert(self, candidate_df: pd.DataFrame):
        """Top level method to make avro alert.

        Args:
            candidate_df (pandas.core.DataFrame): dataframe of candidates.
        """
        t_start = time.time()

        logger.debug(
            f"Avro Packet Maker: parsing {len(candidate_df)} candidates from dataframe"
        )
        # Read dataframe into dictionary of candidates
        all_cands = self.read_input_df(candidate_df)

        num_cands = len(all_cands)  # total candidates to process
        successes = 0  # successful kafka producing of a candidate

        # TODO topic_name should match night of images
        # # Calculate the alert_date from list of jds for this night
        # jd_list = np.array([o['jd'] for o in all_cands])
        # minjd = np.min(jd_list)
        # alert_date = Time(minjd, format = 'jd').tt.datetime.strftime('%Y%m%d')
        # topic_name = f'winter_{alert_date}'
        # logger.info(f'Sending to topic: {topic_name}')

        # TODO topic_name creation used for testing
        topic_name = f"winter_{datetime.datetime.utcnow().strftime('%Y%m%d')}"
        logger.debug(f"topic name {topic_name}")

        self.get_sub_output_dir().mkdir(exist_ok=True, parents=True)

        # Send/make avro packet for each candidate from the dataframe
        for cand in all_cands:
            flag = self._make_ind_packet(cand, self.schema, topic_name)
            if flag > 0:
                successes += 1

            # used for testing script below for checking/opening avro packets
            # cand_id = cand["candid"] # save candid of avro packet to open

        t_end = time.time()
        logger.debug(
            f"Took {(t_end - t_start):.2f} seconds to process {num_cands} candidates. "
            f"{successes} of {num_cands} successfully {self._success_message()}."
        )
