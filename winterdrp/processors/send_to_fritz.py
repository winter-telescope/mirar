import logging
import requests
import pandas as pd
import os, gzip, io

import astropy
from astropy.time import Time
from astropy.io import ascii


from winterdrp.processors.base_processor import BaseDataframeProcessor
from winterdrp.paths import get_output_dir

logger = logging.getLogger(__name__)

class SendToFritz(BaseDataframeProcessor):
    def __init__(self, 
                output_sub_dir: str, 
                token = 'b8d92664-f8ab-407d-b578-648d7237d437',
                group_ids = [1431],
                base_name = 'WIRC',
                *args,
                **kwargs):
        super(SendToFritz, self).__init__(*args, **kwargs)
        self.token = token
        self.group_ids = group_ids
        self.base_name = base_name

    def _apply_to_candidates(
            self,
            candidate_table: pd.DataFrame,
    ) -> pd.DataFrame:
        logger.info("In SendToFritz")
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

    def generate_dict(self, cand):
        """Make dictionary as per API formatting"""
        formatted_dict = {}
        for cand_key in cand.keys():
            if cand_key == 'ra':
                formatted_dict['ra'] = cand[cand_key]
            if cand_key =='dec':
                formatted_dict['dec'] = cand[cand_key]

        formatted_dict['id'] = cand['objectId']
        formatted_dict['group_ids'] = self.group_ids
        logger.info(f"dict:{formatted_dict}")

        return formatted_dict

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

    def api(self, method, endpoint, data=None):
        headers = {'Authorization': f'token {self.token}'}
        response = requests.request(method, endpoint, json=data, headers=headers)
        return response

    def add_new_source(self, cand):
        data = self.generate_dict(cand)
        response = self.api('POST', 'https://fritz.science/api/sources', data)
        return response
    
    def upload_thumbnail(self, cand):
        """Post new thumbnail to Fritz.

        Format of thumbnail payload:
        { "obj_id": "string",  "data": "string",  "ttype": "string"}
        """
        fritz_to_cand = {"new": 'SciBitIm', "ref": 'RefBitIm', "sub": 'DiffBitIm'}
        thumbnail_payload = {}
        thumbnail_payload["obj_id"] = cand['objectId']

        buf = io.BytesIO(cand['SciBitIm']) 
        gzip_f = gzip.GzipFile(fileobj=buf)
        content = gzip_f.read() 
        logger.info(f'content type{type(content)}')       
        thumbnail_payload["data"] = content.decode('latin-1')
        logger.info(f'data type to string: {type(thumbnail_payload["data"])}')
        thumbnail_payload["ttype"] = "new"

        # for fritz_key in fritz_to_cand.keys():
        #     cand_key = fritz_to_cand[fritz_key]
        #     logger.info(f'data: {cand[cand_key]}')
        #     logger.info(f'data type: {type(cand[cand_key])}')
        #     logger.info(f'data type to string: {type(cand[cand_key].)}')


        #     buf = io.BytesIO(cand[cand_key]) 
        #     gzip_f = gzip.GzipFile(fileobj=buf)
        #     content = gzip_f.read()
        #     logger.info(f"data type 2: {type(content)}")

        #     content = content.decode('latin-1')
        #     logger.info(f"data type 3: {type(content)}")

        #     thumbnail_payload['data'] = content
            
        #     thumbnail_payload['ttype'] = fritz_key
        #     logger.info(f'{fritz_key}: {type(thumbnail_payload["data"])}')

        response = self.api('POST', 'https://fritz.science/api/thumbnail', thumbnail_payload)
        logger.info(f'candid {cand["objectId"]}: {thumbnail_payload["ttype"]}, thumbnail response:{response}')

    def make_alert(self, cand_table):
        all_cands = self.read_input_df(cand_table)

        last_name = None
        cand_id = 700
          
    
        for cand in all_cands:
            cand_jd = cand['jd']
            cand_name = self.get_next_name(last_name, str(cand_jd))

            #TODO candid should be coming from naming database
            cand['candid'] = cand_id
            cand_id += 1
            
            # TODO check if cand is new? 
            # cand_name, new_status = self.check_and_insert_source(cand_name, cand)

            #TODO remove once new_status is up-to-date: dummy line
            new_status = True

            if new_status:
                last_name = cand_name    
            cand['objectId'] = cand_name

            
            source_response = self.add_new_source(cand)
            logger.info(source_response)
            thumbnail_response = self.upload_thumbnail(cand)

            break
