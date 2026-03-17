<p align="center">
  <img src="assets/DirectClean_Logo.png" width="250" alt="DirectClean Logo">
</p>

# DirectClean

A comprehensive preprocessing pipeline for **Oxford Nanopore (ONT) direct-cDNA sequencing data**. DirectClean integrates strand orientation, artifact removal, and chimeric read rescue to produce clean, properly oriented FASTQ files optimized for downstream transcript quantification and gene fusion analysis.

**At a glance:**
* **What it removes**: foldback inversion reads (self-inverted artifacts) and reads that cannot be strand-oriented (missing primer signals).
* **What it rescues** (chopped at the artifact junction, flanking sub-reads kept): reads containing internal TSO/RTP adapter junctions (concatemers from ligation) and reads containing homopolymer-mediated RT template switching junctions.

## Motivation

While Pychopper is the standard tool for orienting and rescuing full-length ONT cDNA reads, direct-cDNA library preparation inherently introduces artifacts that Pychopper does not address:

* **Foldback inversions**: The sequenced strand folds back on itself during reverse transcription, producing a self-inverted chimeric read.
* **Homopolymer-mediated RT template switching**: The reverse transcriptase detaches at A/T-rich regions and re-primes on unrelated transcripts. These artifacts mimic genuine gene fusions and inflate false-positive rates in downstream analyses.

DirectClean solves this by combining established tools (Breakinator, Restrander) with novel adapter and homopolymer rescue algorithms into a single, end-to-end pipeline. 

## Pipeline Architecture & Algorithm

Stages 1 and 2 discard definitively artifactual or unorientable reads. Stages 3, 4, and 5 function as rescue modules; they identify internal artifact junctions, chop the chimeric sequences, and rescue the valid flanking sub-reads.

| Stage | Module | Action | Description |
| :--- | :--- | :--- | :--- |
| 1 | Breakinator | Filter | Removes foldback inversion artifacts. |
| 2 | Restrander | Orient & Filter | Orients reads to 5'→3' and removes invalid primer configurations (e.g., RTP-RTP, TSO-TSO). |
| 3 | Unknowns Rescue | Rescue | Recovers orientable reads from Stage 2 unknowns via internal adapter detection. |
| 4 | Adapter Rescue | Rescue | Detects internal TSO/RTP adapters in oriented reads, chops, and rescues sub-reads. |
| 5 | Homopolymer Filter | Rescue | Identifies RT template switching at A/T-rich chimeric junctions. |

Homopolymer Detection Algorithm: Following a splice-aware alignment via minimap2, DirectClean evaluates chimeric junctions defined by supplementary alignments. A 10 bp sliding window scans the sequence flanking each junction. A junction is classified as an RT template switching artifact if any window meets both criteria:
1. A/T base density ≥ 0.85
2. Longest consecutive A or T run ≥ 5 bp

Reads containing flagged junctions are chopped, and sub-reads ≥ 100 bp are rescued.

## Installation

```bash
# 1. Create the conda environment with all dependencies
mamba env create -f environment.yml
mamba activate directclean

# 2. Install DirectClean
poetry install
```

*(Note: External dependencies including minimap2, samtools, breakinator, and restrander are resolved within the conda environment.)*

## Usage

```bash
directclean \
  -i raw_reads.fastq \
  -r genome.fa \
  -o results/ \
  -t 8 \
  -j gencode.v41.bed12
```

*Providing a junction BED file (-j, e.g., GENCODE BED12) is highly recommended for guided splice-aware alignment.*

### Key Parameters

| Argument | Default | Description |
| :--- | :--- | :--- |
| `-i`, `--input` | Required | Raw FASTQ from ONT direct-cDNA sequencing. |
| `-r`, `--reference`| Required | Reference genome FASTA. |
| `-o`, `--output` | Required | Output directory path. |
| `-t`, `--threads` | 4 | Threads allocated for minimap2, samtools, and breakinator. |
| `-j`, `--junc-bed` | None | Junction BED12 file to guide alignment. |
| `--density-threshold`| 0.85 | A/T density threshold for homopolymer detection. |
| `--min-run` | 5 | Minimum consecutive A/T run length for detection. |

*Run `directclean --help` for the complete list of available options.*

## Output Structure

```text
results/
├── directclean.cleaned.fastq          # Primary output: clean reads + rescued sub-reads
├── directclean.rescued.fastq          # Sub-reads specifically rescued via homopolymer chopping
├── directclean.homopolymer_report.tsv # Per-read artifact classification and breakpoints
├── intermediates/                     # FASTQ/BAM files generated after each distinct stage
└── reports/                           # Detailed statistics, including Stage 4 adapter rescue metrics
```

The primary downstream file is `directclean.cleaned.fastq`, which can be directly used for transcript quantification (e.g., IsoQuant, FLAIR) and gene fusion calling (e.g., FusionSeeker, JAFFAL).

## Performance

Benchmarked on 5.35M direct-cDNA reads from the VCaP prostate cancer cell line:

| Metric | Pychopper | DirectClean |
| :--- | :--- | :--- |
| Retention rate | 57.6% | 65.3% |
| FSM isoforms detected | 17,873 | 20,535 |
| Validated fusions detected (n=99)| 37 | 49 |
| Residual homopolymer artifacts | 70,140 | 0 |

## Citation
Manuscript in preparation.

## License
MIT

## Contact
* Qingxiang Guo - qingxiang.guo@northwestern.edu
* Rendong Yang Lab - https://github.com/ylab-hi