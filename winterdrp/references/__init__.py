from winterdrp.errors import ProcessorError
from winterdrp.references.base_reference_generator import BaseReferenceGenerator
from winterdrp.references.ps1 import PS1Ref
from winterdrp.references.sdss import SDSSRef
from winterdrp.references.wirc import WIRCRef


class ReferenceImageError(ProcessorError):
    pass
