# DirectClean

A comprehensive pipeline for removing RT artifacts from Oxford Nanopore Direct-cDNA sequencing data.

## Overview

DirectClean integrates Breakinator, Restrander, and custom filtering algorithms to clean Direct-cDNA reads through a complete workflow:

- **Breakinator** removes foldback inversion artifacts
- **Restrander** corrects strand orientation and removes RTP-RTP/TSO-TSO artifacts  
- **Adapter rescue** detects and removes internal TSO/RTP sequences
- **Homopolymer filtering** removes template switching artifacts at poly-A/T junctions

The cleaned output is optimized for fusion detection and multi-gene transcript analysis.

## Installation

## Usage
```bash
directclean \
    --input raw_reads.fastq \
    --reference genome.fa \
    --output results/ \
    --threads 8
```

## Output

- `cleaned.fastq` - Filtered reads ready for analysis
- `removed.fastq` - Artifact reads
- `report.tsv` - Filtering statistics

## Citation

Manuscript in preparation.

## License

MIT License
