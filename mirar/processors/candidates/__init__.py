"""
Central module for candidate detection and extraction.
"""
from mirar.processors.candidates.candidate_detector import DetectCandidates
from mirar.processors.candidates.candidate_extractor import (
    SourceTablefromCoordinates,
    SourceTablefromHeader,
)
from mirar.processors.candidates.candidate_filter import FilterCandidates
from mirar.processors.candidates.dataframe_writer import DataframeWriter
from mirar.processors.candidates.edge_mask import EdgeCandidatesMask
from mirar.processors.candidates.namer import CandidateNamer
