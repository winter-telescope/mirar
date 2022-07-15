import logging
from multiprocessing import get_context

import astropy.table
from astropy.time import Time
from astropy.io import ascii
import pandas as pd

import io, gzip, os, sys
import avro, fastavro
from avro import schema
import avro.schema
import avro.io
import avro.tool
from avro.datafile import DataFileWriter, DataFileReader
from avro.io import DatumWriter, DatumReader
import confluent_kafka
import copy
import json
import time

from winterdrp.processors.base_processor import BaseDataframeProcessor
from winterdrp.paths import get_output_dir

logger = logging.getLogger(__name__)


class AvroPacketMaker(BaseDataframeProcessor):
    """Class to generate Avro Packets from a dataframe of candidates.

    Attributes:
        output_sub_dir (str): output data path.
        base_name (str): 4-letter code for telescope.
    """

    def __init__(self, 
                output_sub_dir: str, 
                base_name: str,
                save_local=False,
                use_database=False,
                broadcast=True,
                *args,
                **kwargs):
        super(AvroPacketMaker, self).__init__(*args, **kwargs)
        self.output_sub_dir = output_sub_dir
        self.base_name = base_name
        self.save_local = save_local
        self.use_database = use_database
        self.broadcast = broadcast

    def _apply_to_candidates(
            self,
            candidate_table: pd.DataFrame,
    ) -> pd.DataFrame:
        self.make_alert(candidate_table) 
        return candidate_table

    def read_input_df(self, df):
        """Takes a DataFrame, which has multiple candidate 
        and creates list of dictionaries, each dictionary 
        representing a single candidate.

        Args:
            df (pandas.core.frame.DataFrame): dataframe of all candidates.
        
        Returns:
            (list[dict]): list of dictionaries, each a candidate.
        """
        all_candidates = []   
        
        for i in range(0, len(df)):
            candidate = {} 
            for key in df.keys():
                try: 
                    if type(df.iloc[i].get(key)) is str:
                        candidate[key] = df.iloc[i].get(key)   
                    else:
                        # change to native python type
                        candidate[key] = df.iloc[i].get(key).item()
                except AttributeError: # for IOBytes objs
                    candidate[key] = df.iloc[i].get(key).getvalue()
                                                 
            all_candidates.append(candidate)

        return all_candidates 

    def combine_schemas(self, schema_files):
        """Combine multiple nested schemas into a single schema.
        Modified from https://github.com/dekishalay/pgirdps

        Args:
            schema_files (list[str]): list of file paths to .avsc schemas.

        Returns:
            (dict): built avro schema.
        """
        known_schemas = avro.schema.Names() # avro.schema.Names object
        
        for s in schema_files:
            schema = self.load_single_avsc(s, known_schemas) 

        # using schema.to_json() doesn't fully propagate the nested schemas
        # work around as below
        props = dict(schema.props)
        fields_json = [field.to_json() for field in props['fields']]
        props['fields'] = fields_json

        return props

    def load_single_avsc(self, file_path, names):
        """Load a single avsc file.
        Modified from https://github.com/dekishalay/pgirdps

        Args:
            file_path (str): file path of the .avsc schema to build.
            names (avro.schema.Names): an avro schema.

        Returns:
            (avro.schema.RecordSchema): data in avro schema format.
        """
        curdir = os.path.dirname(__file__)
        file_path = os.path.join(curdir, file_path)

        with open(file_path) as file_text:
            json_data = json.load(file_text)

        # SchemaFromJSONData not working
        # ##### works only with avro version 1.10.1 #####
        schema = avro.schema.SchemaFromJSONData(json_data, names) 

        return schema

    def get_next_name(self, lastname, candjd, bwfile = 'badwords.txt', begcount = 'aaaaaaa'):
        """Creates candidate name following the naming format of 'WNTR22aaaaaaa' .
        Modified from https://github.com/dekishalay/pgirdps

        Args:
            lastname(str): last used candidate name.
            candjd (str): candidate's JD.
            bwfile (str): file name of .txt file of excluded words. 
            begcount (str): string to start naming convention.
        
        Returns:
            (str): next candidate name in sequence.
        """ 
        curyear = Time(candjd, format = 'jd').datetime.strftime('%Y')[2:4]

        if lastname is None:
            #If this is the first source being named
            newname = self.base_name + curyear + begcount
            return newname
        
        lastyear = lastname[4:6]

        if curyear != lastyear:
            #If this is the first candidate of the new year, start with aaaaa
            newname = self.base_name + curyear + begcount
            return newname
        else:
            lastcount = lastname[6:]
            charpos = len(lastcount) - 1
            # will iteratively try to increment characters starting from the last
            inctrue = False
            usestring = ''
            while charpos >= 0:
                cref = lastcount[charpos]
                if inctrue:
                    usestring = cref + usestring
                    charpos -= 1
                    continue
                creford = ord(cref)
                #increment each character, if at 'z', increment the next one
                if creford + 1 > 122:
                    usestring = 'a' + usestring
                else:
                    nextchar = chr(creford+1)
                    usestring = nextchar + usestring
                    inctrue = True
                charpos -= 1
                continue
            
            newname = self.base_name + curyear + usestring

            curdir = os.path.dirname(__file__) # /data/sulekha/winterdrp/winterdrp/processors/alert_packets
            file_path = os.path.join(curdir, bwfile)

            bwlist = ascii.read(file_path, format = 'no_header')
            isbw = False
            # check for bad word
            for i in range(len(bwlist)):
                if usestring.find(str(bwlist['col1'][i])) != -1:
                    #Shame shame
                    isbw = True
                    break
            if isbw:
                # increment the name with a recursive call
                return self.get_next_name(newname, candjd)
            else:
                return newname

    def write_avro_data(self, json, avro_schema):
        """Encode json into avro format given a schema.

        Args:
            json (dict): data to put into schema, can be dict.
            avro_schema (avro.schema.RecordSchema): output schema.

        Returns:
            (io.BytesIO object): input json in schema format.
        """
        writer = avro.io.DatumWriter(avro_schema)
        bytes_io = io.BytesIO()
        encoder = avro.io.BinaryEncoder(bytes_io)
        writer.write(json, encoder)
        return bytes_io

    def read_avro_data(self, bytes_io, avro_schema):
        """Read avro data and decode with a given schema.
        
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

    def _send_alert(self, topicname, records, schema):
        """Send an avro "packet" to a particular topic at IPAC.
        Modified from: https://github.com/dekishalay/pgirdps

        Args:
            topicname (str): name of the topic sending to, e.g. ztf_20191221_programid2_zuds.
            records (list): a list of dictionaries (the avro packet to send).
            schema (dict): schema definition.
        """
        out = io.BytesIO()
        
        fastavro.writer(out, schema, records)
        out.seek(0) # go back to the beginning
        
        # Connect to the IPAC Kafka brokers
        producer = confluent_kafka.Producer({'bootstrap.servers': 'ztfalerts04.ipac.caltech.edu:9092,ztfalerts05.ipac.caltech.edu:9092,ztfalerts06.ipac.caltech.edu:9092'})

        # Send an avro alert
        producer.produce(topic=topicname, value=out.read())
        producer.flush()

    def create_alert_packet(self, cand, scicut, refcut, diffcut, cm_radius=8.0, search_history=90.0):
        """Create top level avro packet from input candidate.
        Args:
            cand (dict): all data of a single candidate.
            scicut (bytes): science image cutout.
            refcut (bytes): reference image cutout.
            diffcut (bytes): difference image cutout.
            cm_radius (float): cross-match radius in arcsec.
            search_history (float): in days.

        Returns:
            (dict): schema in dictionary format.
        """
        # TODO populate the candidate history (prv_candidate)
        prev_cands = []
        # logger.info(f' {len(prev_cands)} prev candidates found ########')

        alert = {"schemavsn": "0.1", "publisher": "winter_test", 
		"cutoutScience": scicut,
		"cutoutTemplate": refcut,
		"cutoutDifference": diffcut,
		"objectId": cand['objectId'],
		"candid": cand['candid'], 
		"candidate": cand,
		"prv_candidates": prev_cands
		}

        return alert
    
    def get_sub_output_dir(self):
        """Returns path of output subdirectory."""
        return get_output_dir(self.output_sub_dir, self.night_sub_dir)

    def _make_avro_output_dir(self):
        """ Make 'avro' subdirectory for avro packet output. """
        avro_output_dir = self.get_sub_output_dir()
        try: # subdir doesn't exist
            os.makedirs(avro_output_dir, exist_ok=True)
        except OSError:
            pass


    def _save_local(self, candid, records, schema):
        """Save avro file in given schema to output subdirectory.
        
        Args:
            candid (number type): candidate id.
            records (list): a list of dictionaries
            schema (avro.schema.RecordSchema): schema definition.
        """
        self._make_avro_output_dir()        
        avro_packet_path = os.path.join(self.get_sub_output_dir(), str(candid) + '.avro')
        out = open(avro_packet_path, 'wb')
        # logger.info(f'out file: {out}')
        fastavro.writer(out, schema, records)
        out.close()

    def save_alert_packet(self, packet, cand, schema):
        """Saves packet as .avro to output subdirectory.
        
        Args:
            packet (dict): candidate data in avro packed dict format.
            cand (dict): all data of a single candidate.
            schema (avro.schema.RecordSchema): schema definition.

        Returns:
            (int): 1 if save successful. Else, -1.
        """        
        try:
            self._save_local(cand['candid'], [packet], schema)
            # logger.info(f"Saved candid {cand['candid']}: {cand['objectId']}")
            return 1
        except Exception as e:
            logger.info(f'{e}')
            logger.info(f"Could not save candid {cand['candid']}")
            return -1
    
    def broadcast_alert_packet(self, packet, cand, schema, topic_name, cand_num, num_cands):
        """
        Sends avro-formatted packets to specified topicname using Kafka. 
        Modified from https://github.com/dekishalay/pgirdps

        Args:
            packet (dict): candidate data in avro packed dict format.
            cand (dict): all data of a single candidate.
            schema (dict): schema definition.
            topic_name (str): name of the topic sending to, e.g. ztf_20191221_programid2_zuds.
            cand_num (int): number of current candidate being sent.
            num_cands (int): total number of candidates to send. 

        Returns:
            (int): 1 if broadcast successful or -1 if candidate not sent.
        """
        try:
            self._send_alert(topic_name, [packet], schema)
            # logger.info(f"Sent candid {cand['candid']}, name {cand['objectId']}, {cand_num} out of {num_cands}")
            return 1
        except OSError:
            logger.info(f"Could not send candid {cand['candid']}")
            return -1

    def _make_ind_packet(self, cand, schema, topic_name, cand_num, num_cands):
        """Makes a single avro alert. Saves or broadcasts alert based on global var.
        
        Args:
            cand (dict): all data of a single candidate.
            schema (dict): schema definition.
            topic_name (str): name of the topic sending to, e.g. ztf_20191221_programid2_zuds.
            cand_num (int): number of current candidate being sent.
            num_cands (int): total number of candidates to send. 

        Returns:
            (int): 1 if broadcast successful,
                   0 if alert neither saved nor broadcast,
                  -1 if candidate not sent or saved.
        """
        # Cutouts are include in the top level alert schema
        scicut = cand.pop('SciBitIm')
        refcut = cand.pop('RefBitIm')
        diffcut = cand.pop('DiffBitIm')

        packet = self.create_alert_packet(cand, scicut, refcut, diffcut)

        flag = 0
        if self.save_local:
            flag += self.save_alert_packet(packet, cand, schema)
        if self.broadcast:
            flag += self.broadcast_alert_packet(packet, cand, schema, 
                                                    topic_name, cand_num, num_cands)
        
        return flag


    def make_alert(self, df=None):
        """Top level method to make avro alert.
        
        Args:
            df (pandas.core.DataFrame): dataframe of candidates.
        """
        t0 = time.time()       
        
        if not self.use_database and df is not None:
            # input dataframe to avro_creation processor
            logger.info(f'Avro Packet Maker: parsing {len(df)} candidates from dataframe')
            all_cands = self.read_input_df(df)
        
        num_cands = len(all_cands)
        successes = 0

        # Make top level alert schema
        schema = self.combine_schemas(["alert_schema/candidate.avsc", 
                                    "alert_schema/prv_candidate.avsc", 
                                    "alert_schema/alert.avsc"])

        # TODO: fake!! remove; lastname, cand_id needs to updated from database
        last_name = None
        cand_id = 700
        cand_num = 1 # for keeping track of broadcast alerts

        for cand in all_cands:
            cand_jd = cand['jd']
            cand_name = self.get_next_name(last_name, str(cand_jd))

             # TODO create alert_date once (before loop)from list of jd_list 
            # Calculate the alert_date for this night
            # alert_date = Time(cand_jd, format = 'jd').tt.datetime.strftime('%Y%m%d')
            # logger.info(alert_date)
            # topic_name = 'winter_%s'%alert_date
            topic_name = 'winter_20220714'

            # TODO candid should be coming from naming database
            cand['candid'] = cand_id
            cand_id += 1
            
            # TODO check if cand is new? 
            # cand_name, new_status = self.check_and_insert_source(cand_name, cand)

            #TODO remove once new_status is up-to-date: dummy line
            new_status = True

            if new_status:
                last_name = cand_name    
            cand['objectId'] = cand_name
            
            flag = self._make_ind_packet(cand, schema, topic_name, cand_num, num_cands)
            cand_num += 1
            if flag > 0:
                successes += 1
           
        t1 = time.time()
        logger.info('###########################################')
        logger.info(f"Took {(t1 - t0):.2f} seconds to process {num_cands} candidates")
        logger.info(f'{successes} of {num_cands} successfully sent')
        logger.info('###########################################')


        # Read data from an avro file
        last_file = str(cand_id - 1) + '.avro'
        last_packet_path = os.path.join(self.get_sub_output_dir(), last_file)
        with open(last_packet_path, 'rb') as f:
            reader = DataFileReader(f, DatumReader())
            metadata = copy.deepcopy(reader.meta)
            schema_from_file = json.loads(metadata['avro.schema'])
            cand_data = [field_data for field_data in reader]
            reader.close()
        
        ## cand_data is a list with the first item as the nested
        ## dictionary of the avro schema
        ## length of cand_data: 1
        ## 0-index item: dictionary, length 9, keys:
        ## ['schemavsn', 'publisher', 'objectId', 'candid', 'candidate', 
        ##  'prv_candidates', 'cutoutScience', 'cutoutTemplate', 'cutoutDifference'])]

        cand_dict = cand_data[0]
        single_cand = cand_dict['candidate']
        jd_val = single_cand['jd']
        logger.info(f'jd: {jd_val}')

        # logger.info(f'Schema that we parsed:\n {schema}')
        # logger.info(f'Schema from candidate .avro file:\n {schema_from_file}')
        # logger.info(f'Candidate:\n {cand_data}')
        # logger.info(f'type{type(cand_data)}')