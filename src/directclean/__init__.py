"""
DirectClean — Remove RT artifacts from Oxford Nanopore Direct-cDNA sequencing.

A comprehensive preprocessing pipeline that detects and removes:
  1. Internal TSO/RTP adapter chimeras (Rescuer module)
  2. Homopolymer-mediated RT template switching artifacts (Filter module)

Typical CLI usage::

    directclean -i reads.fastq -r genome.fa -o results/ -t 8
"""

__version__ = "0.1.0"
__author__ = "Qingxiang Guo"