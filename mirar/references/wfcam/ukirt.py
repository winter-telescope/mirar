import logging

from astroquery.ukidss import UkidssClass
from astroquery.wfau import BaseWFAUClass
from astrosurveyutils.surveys import MOCSurvey

from mirar.references.wfcam.utils import find_ukirt_surveys
from mirar.references.wfcam.wfcam import WFCAMStackedRefOnline

logger = logging.getLogger(__name__)


class UKIRTStackRefOnline(WFCAMStackedRefOnline):
    def get_surveys(self, ra: float, dec: float) -> list[MOCSurvey]:
        return find_ukirt_surveys(ra=ra, dec=dec, band=self.filter_name)

    def get_query_class(self) -> BaseWFAUClass:
        return UkidssClass()
