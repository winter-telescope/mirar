"""
Central module for candidate detection and extraction.
"""
from mirar.processors.sources.candidate_filter import BaseSourceFilter
from mirar.processors.sources.dataframe_writer import DataframeWriter
from mirar.processors.sources.edge_mask import EdgeSourcesMask
from mirar.processors.sources.namer import CandidateNamer
from mirar.processors.sources.source_detector import SourceDetector
from mirar.processors.sources.source_table_builder import (
    SourceTablefromCoordinates,
    SourceTablefromHeader,
)
