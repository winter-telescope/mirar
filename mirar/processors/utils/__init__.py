"""
Module for general utility processors such as I/O and interacting with metadata
"""
from mirar.processors.utils.header_annotate import HeaderAnnotator
from mirar.processors.utils.header_reader import HeaderReader
from mirar.processors.utils.image_loader import ImageLoader, MEFLoader
from mirar.processors.utils.image_modifier import CustomImageModifier
from mirar.processors.utils.image_saver import ImageSaver
from mirar.processors.utils.image_selector import (
    ImageBatcher,
    ImageDebatcher,
    ImageSelector,
    select_from_images,
)
from mirar.processors.utils.multi_ext_parser import MultiExtParser
