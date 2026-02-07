"""
Homopolymer-mediated RT artifact filter.

This sub-package detects and removes chimeric reads caused by
reverse-transcriptase template switching at poly-A/T regions
in Oxford Nanopore Direct-cDNA sequencing data.

Public API::

    from directclean.filter import ArtifactClassifier, HomopolymerConfig

    classifier = ArtifactClassifier(
        bam_path="aligned.bam",
        input_fastq="reads.fastq",
        output_dir="output/",
        config=HomopolymerConfig(scan_window=10, density_threshold=0.8, min_run=3),
    )
    report = classifier.run()
"""

from directclean.filter.artifact_classifier import ArtifactClassifier, FilterReport
from directclean.filter.homopolymer import (
    HomopolymerConfig,
    HomopolymerDetector,
    JunctionVerdict,
    ReadVerdict,
)
from directclean.filter.junction_parser import (
    ChimericRead,
    JunctionInfo,
    SegmentInfo,
    parse_chimeric_read,
    iter_chimeric_reads,
)

__all__ = [
    "ArtifactClassifier",
    "FilterReport",
    "HomopolymerConfig",
    "HomopolymerDetector",
    "JunctionVerdict",
    "ReadVerdict",
    "ChimericRead",
    "JunctionInfo",
    "SegmentInfo",
    "parse_chimeric_read",
    "iter_chimeric_reads",
]