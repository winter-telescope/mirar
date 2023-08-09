"""
Central module for candidate detection and extraction.
"""
from mirar.processors.sources.edge_mask import EdgeSourcesMask
from mirar.processors.sources.namer import CandidateNamer
from mirar.processors.sources.source_detector import ZOGYSourceDetector
from mirar.processors.sources.source_exporter import SourceWriter
from mirar.processors.sources.source_filter import BaseSourceFilter
from mirar.processors.sources.source_loader import SourceLoader
from mirar.processors.sources.source_table_builder import (
    SourceTablefromCoordinates,
    SourceTablefromHeader,
)
from mirar.processors.sources.source_table_modifier import CustomSourceModifier
