# DirectClean

**Artifact-aware preprocessing for Oxford Nanopore direct-cDNA
sequencing**

DirectClean is a preprocessing toolkit for ONT direct-cDNA sequencing
reads. It integrates strand correction, artifact detection, and sequence
rescue to improve recovery of usable transcript sequences for transcript
discovery, isoform analysis, and fusion detection.

Unlike conventional read filtering approaches that discard
artifact-containing reads, DirectClean aims to preserve informative
sequence by identifying artifact structures and applying targeted
trimming, splitting, and rescue strategies.

------------------------------------------------------------------------

## Overview

Oxford Nanopore direct-cDNA sequencing provides long transcript
molecules but can contain multiple artifact classes introduced during
library preparation and reverse transcription.

DirectClean addresses these challenges through a multi-stage workflow:

``` mermaid
flowchart TD
    A[Raw ONT direct-cDNA FASTQ] --> B[Stage 1<br/>Breakinator<br/>Foldback artifact detection]
    B --> C[Stage 2<br/>Restrander<br/>Strand correction]
    C --> D[Stage 3<br/>Unknowns Rescue<br/>Recover usable reads]
    D --> E[Stage 4<br/>Adapter Rescue<br/>Remove adapters and split supported concatemers]
    E --> F[Stage 5<br/>Homopolymer Rescue<br/>Resolve RT template-switch artifacts]
    F --> G[Clean FASTQ<br/>Transcript and fusion analysis]
```

  -----------------------------------------------------------------------
  Stage                   Component               Description
  ----------------------- ----------------------- -----------------------
  1                       Breakinator             Detects and removes
                                                  foldback inversion
                                                  artifacts

  2                       Restrander              Corrects read
                                                  orientation and
                                                  identifies
                                                  strand-related
                                                  artifacts

  3                       Unknowns Rescue         Recovers usable
                                                  sequence from reads
                                                  with unresolved
                                                  classification

  4                       Adapter Rescue          Removes
                                                  adapter-associated
                                                  artifacts and splits
                                                  supported concatemer
                                                  structures

  5                       Homopolymer Rescue      Detects and resolves
                                                  homopolymer-mediated RT
                                                  template-switch
                                                  artifacts
  -----------------------------------------------------------------------

------------------------------------------------------------------------

# Installation

## Recommended: Bioconda

The recommended installation method is through Bioconda:

``` bash
mamba create -n directclean \
    -c conda-forge \
    -c bioconda \
    directclean

mamba activate directclean
```

Verify installation:

``` bash
directclean --help
```

DirectClean automatically installs required dependencies, including:

-   minimap2
-   samtools
-   breakinator
-   restrander

------------------------------------------------------------------------

# Quick usage

Example:

``` bash
directclean \
    -i raw_reads.fastq \
    -r genome.fa \
    -o results \
    -t 8 \
    -j annotation.bed12
```

Main output:

    results/
    ├── directclean.cleaned.fastq
    ├── directclean.report.html
    └── intermediates/

The cleaned FASTQ can be directly used for downstream transcriptome
analysis.

------------------------------------------------------------------------

# Why DirectClean?

[Pychopper](https://github.com/epi2me-labs/pychopper) is widely used for
ONT direct-cDNA preprocessing, particularly for primer detection and
strand orientation correction.

DirectClean complements this functionality by addressing additional
artifact classes generated during direct-cDNA sequencing and by rescuing
usable sequence from artifact-containing reads.

## Artifact classes addressed

  Capability                                                  Pychopper   DirectClean
  ---------------------------------------------------------- ----------- -------------
  Primer-based orientation                                        ✓            ✓
  Strand correction                                               ✓            ✓
  Foldback inversion artifact detection                          --            ✓
  Adapter-associated concatemer resolution                       --            ✓
  Homopolymer-mediated RT template-switch detection              --            ✓
  Rescue of usable sequence from artifact-containing reads       --            ✓

------------------------------------------------------------------------

# Benchmark against Pychopper

DirectClean was evaluated on an ONT direct-cDNA dataset from VCaP
prostate cancer cells.

  Metric                  Pychopper   DirectClean
  --------------------- ----------- -------------
  Input reads             5,348,910     5,348,910
  Output-record yield         57.6%     **65.3%**

DirectClean recovered more output records while performing additional
artifact-resolution steps.

Importantly, output-record yield alone does not define transcript
retention. DirectClean performs artifact-aware trimming and splitting to
preserve usable transcript sequence rather than simply retaining or
discarding complete reads.

------------------------------------------------------------------------

# Artifact resolution

## Foldback artifacts

Foldback structures can arise when a molecule contains inverted sequence
copies caused by strand folding events.

DirectClean uses Breakinator-based detection to identify and remove
these artifacts before downstream analysis.

## Adapter-associated structures

Direct-cDNA reads may contain adapter-derived structures or
concatemer-like artifacts.

DirectClean distinguishes these cases by:

-   removing unsupported adapter sequence;
-   splitting supported internal concatemer junctions;
-   retaining resulting transcript fragments when sequence evidence
    supports recovery.

## Homopolymer-mediated RT template switching

Reverse transcription template switching can generate artificial
junctions associated with A/T-rich homopolymer regions.

DirectClean identifies these events using sequence-context features and
rescues supported subreads.

------------------------------------------------------------------------

# HTML report

DirectClean generates a self-contained interactive HTML report
containing:

-   executive summary statistics;
-   read-flow waterfall visualization;
-   stage-specific classification summaries;
-   artifact-resolution statistics.

This allows users to inspect processing outcomes without additional
visualization software.

------------------------------------------------------------------------

# Output files

Typical output structure:

    results/
    ├── directclean.cleaned.fastq
    ├── directclean.report.html
    ├── statistics/
    └── intermediates/

`directclean.cleaned.fastq` is recommended for:

-   transcript reconstruction;
-   isoform quantification;
-   fusion detection;
-   downstream long-read RNA analysis.

------------------------------------------------------------------------

# Developer installation

For development or source installation:

``` bash
git clone https://github.com/ylab-hi/DirectClean.git
cd DirectClean

mamba env create -f environment.yml
mamba activate directclean

pip install -e .
```

------------------------------------------------------------------------

# Citation

If you use DirectClean, please cite:

-   DirectClean manuscript
-   Breakinator
-   Restrander

------------------------------------------------------------------------

# License

MIT
