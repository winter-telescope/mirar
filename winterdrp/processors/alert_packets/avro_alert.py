import logging

import astropy.table
import pandas as pd

import io, gzip, os
import avro, fastavro
from avro import schema
import json
import time

from winterdrp.processors.base_processor import BaseDataframeProcessor
from winterdrp.paths import get_output_dir

logger = logging.getLogger(__name__)


class AvroPacketMaker(BaseDataframeProcessor):
    """Class to generate Avro Packets from a dataframe of candidates.
    Attributes:
        schema_path (str): path for folder with alert.avsc, 
                        candidate.avsc, and prv_candidate.avsc.
        store_data_path (str): path to store avro files in.
    """

    base_key = "avro"

    def __init__(self, 
                output_sub_dir: str = "avro_packets", 
                *args,
                **kwargs):
        super(AvroPacketMaker, self).__init__(*args, **kwargs)
        self.output_sub_dir = output_sub_dir

    def _apply_to_candidates(
            self,
            candidate_table: pd.DataFrame,
    ) -> pd.DataFrame:
        logger.info('in AvroPacketMaker: _apply_to_image')
        self.make_alert(False, candidate_table[0]) # TODO go thru entire list
        return candidate_table
    
    def read_input_df(
        self, 
        df
    ) -> list[dict]:
        """Takes a DataFrame (single table from input list[pd.DataFrame])
        which has multiple candidate & creates list of dictionaries, 
        each dictionary representing a single candidate.

        Args:
            df (pandas.core.frame.DataFrame): dataframe of all candidates.
        
        returns:
            list[dict]: list of dictionaries, each a candidate.
        """
        all_candidates = []   
        
        for i in range(0, len(df)):
            candidate = {} 
            for key in df.keys():
                try: 
                    candidate[key] = df.iloc[i].get(key).item() # change to native python type
                except AttributeError: # for IOBytes objs
                    candidate[key] = df.iloc[i].get(key).getvalue()
            all_candidates.append(candidate)

        return all_candidates 

    def combine_schemas(self, schema_files):
        """Combine multiple nested schemas into a single schema.
        Modified from https://github.com/dekishalay/pgirdps
        """
        known_schemas = avro.schema.Names()
        logger.info(f'known schemas: {known_schemas}')
        #print(known_schemas)
        for s in schema_files:
            schema = self.load_single_avsc(s, known_schemas) # not working 
            #     schema = avro.schema.SchemaFromJSONData(json_data, names)
            # AttributeError: module 'avro.schema' has no attribute 'SchemaFromJSONData'
        # using schema.to_json() doesn't fully propagate the nested schemas
        # work around as below
        logger.info(f'props: {schema.props}')
        logger.info(f'props type: {type(schema.props)}')

        props = dict(schema.props)
        fields_json = [field.to_json() for field in props['fields']]
        props['fields'] = fields_json
        return props

    def load_single_avsc(self, file_path, names):
        """Load a single avsc file.
        Modified from https://github.com/dekishalay/pgirdps
        """
        curdir = os.path.dirname(__file__)
        logger.info(f'cur_dir: {curdir}')
        file_path = os.path.join(curdir, file_path)
        logger.info(f'fil path: {file_path}')

        with open(file_path) as file_text:
            json_data = json.load(file_text)
            logger.info(f'json_data type: {type(json_data)}')
            logger.info(f'json_data dumps type: {type(json.dumps(json_data))}')

        # SchemaFromJSONData not working
        # schema = avro.schema.SchemaFromJSONData(json_data, names)
        schema = avro.schema.parse(json.dumps(json_data))
        logger.info(f'schema type: {type(schema)}')

        return schema

    def write_avro_data(json, avro_schema):
        """Encode json into avro format given a schema.
        Args:
            json (dict): data to put into schema, can be dict.
            avro_schema (avro.schema.RecordSchema): output schema.
        Returns:
            _io.BytesIO object: input json in schema format.
        """
        writer = avro.io.DatumWriter(avro_schema)
        bytes_io = io.BytesIO()
        encoder = avro.io.BinaryEncoder(bytes_io)
        writer.write(json, encoder)
        return bytes_io

    def read_avro_data(bytes_io, avro_schema):
        """Read avro data and decode with a given schema.
        
        Args:
            bytes_io (_io.BytesIO object): data to put into schema, can be dict.
            avro_schema (avro.schema.RecordSchema): schema to decode from.
        Returns:
            dict: schema in dictionary format
        """
        raw_bytes = bytes_io.getvalue()
        bytes_reader = io.BytesIO(raw_bytes)
        decoder = avro.io.BinaryDecoder(bytes_reader)
        reader = avro.io.DatumReader(avro_schema)
        message = reader.read(decoder)
        return message

    def create_alert_packet(self, cand):
        """Create top level avro packet from input candidate.
        Args:
            cand (dict): all data of a single candidate.
        Returns:
            dict: schema in dictionary format.
        """
        # Make candidate subschema

        # TODO populate the candidate history (prv_candidate)
        prev_cands = {}
        logger.info(f'{len(prev_cands)} prev candidates found.')

        # TODO read scicut, refcut, diffcut, candidate, precands in
        # #converting memoryview cutouts to bytes before sending to parallelized broadcast
        # scicut = io.BytesIO(canddict.pop('sci_image')).read()
        # refcut = io.BytesIO(canddict.pop('ref_image')).read()
        # diffcut = io.BytesIO(canddict.pop('diff_image')).read()
        # avro_data = self.write_avro_data(cand, schema)
        # print(f'output: {avro_data}')
        # collected_stamp = self.read_avro_data(avro_data, schema_parsed)

        # TODO fix!
        scicut = cand['SciBitIm'] # should be bytes
        refcut = cand['RefBitIm']
        diffcut = cand['DiffBitIm']

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

    def _save_local(self, candid, records, schema):
        """Save avro file in given schema to output subdirectory.
        
        Args:
            candid (number type): candidate id.
            records (list): a list of dictionaries
            avro_schema (avro.schema.RecordSchema): schema definition.
        """
        avro_sub_dir = '%d.avro'%candid
        filename = get_output_dir(avro_sub_dir, self.output_sub_dir) # possible wrong order
        out = open(filename, 'wb')
        logger.info(f'avro alert saved to {out}')
        fastavro.writer(out, schema, records)
        out.close()

    def save_alert_packet(self, cand, scicut, refcut, diffcut, schema, make_file=False):
        packet = self.create_alert_packet(cand, scicut, refcut, diffcut)
        
        if make_file:
            ## Manual way, equiv to above funcs
            # ## Write data to an avro file
            # filename = self.store_data_path + packet['objectID'] + '.avro'
            # with open(filename, 'wb') as f:
            #     writer = avro.io.DataFileWriter(f, avro.io.DatumWriter(), schema)
            #     writer.append(cand)
            #     writer.close()  \
            try:
                self._save_local(cand['candid'], [packet], schema)
                logger.info('Saved candid %d name %s'%(cand['candid'], cand['objectId']))
            except:
                logger.info('Could not save candid %d'%cand['candid'])
        

    def make_alert(self, useDataBase=False, df=None):
        """Top level method to make avro alert
        Args:
            bytes_io (_io.BytesIO object): data to put into schema, can be dict.
            useDataBase (bool): schema to decode from.
        Returns:
            dict: schema in dictionary format
        """
        t0 = time.time()       
        
        if not useDataBase and df is not None:
            # input dataframe to avro_creation processor
            logger.info(f'{len(df)} candidates in df')
            all_cands = self.read_input_df(df)
            logger.info('####################################')
            logger.info(f'{len(all_cands)} candidates in dict')
        
        ####### Commenting out for push 
        # # TODO change to self.schema_path & string concat
        # schema = self.combine_schemas(["alert_schema/candidate.avsc", 
        #                             "alert_schema/prv_candidate.avsc", 
        #                             "alert_schema/alert.avsc"])

        # # input dataframe to avro_creation processor
        # # df = pd.DataFrame() #TODO fix!
        # # all_cands = self.read_input_df(df)



        num_cands = len(all_cands)
        # logger.info(f'{num_cands} candidates found...making packets')

        # for cand in all_cands:
        #     scicut = cand.pop('SciBitIm')
        #     refcut = cand.pop('RefBitIm')
        #     diffcut = cand.pop('DiffBitIm')
        #     self.save_alert_packet(cand, scicut, refcut, diffcut, schema, True)

        t1 = time.time()
        print('Took %.2f seconds to process %d candidates'%(t1 - t0, num_cands))