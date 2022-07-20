from winterdrp.processors.base_processor import BaseDataframeProcessor
import numpy as np


class CandidateNamer(BaseDataframeProcessor):
    def __init__(self,
                 xmatch_radius_arcsec,
                 *args,
                 **kwargs):
        super(CandidateNamer, self).__init__(*args, **kwargs)
        self.xmatch_radius_arcsec = xmatch_radius_arcsec
        pass
