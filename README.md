<p align="center">
  <img src="assets/DirectClean_Logo.png" width="250" alt="DirectClean Logo">
</p>

# DirectClean

Strand orientation, artifact removal, and chimeric read rescue for Oxford Nanopore direct-cDNA sequencing.

DirectClean processes raw ONT direct-cDNA FASTQ files and produces clean, oriented reads ready for transcript quantification and gene fusion analysis.

**What it removes:** foldback inversion reads (self-inverted artifacts) and reads that cannot be strand-oriented (missing primer signals).

**What it rescues** (chopped at the artifact junction, flanking sub-reads kept): reads containing internal TSO/RTP adapter junctions (concatemers from ligation) and reads containing homopolymer-mediated RT template switching junctions.

## Performance on VCaP direct-cDNA data

Tested on 5.35M reads from the VCaP prostate cancer cell line:

| Metric | Pychopper | DirectClean |
| :--- | :--- | :--- |
| Retention rate | 57.6% | **65.3%** |
| FSM isoforms detected | 17,873 | **20,535** |
| Validated fusions detected (of 99) | 37 | **49** |
| Residual homopolymer artifacts | 70,140 | **0** |

## Why DirectClean?

Oxford Nanopore's [Pychopper](https://github.com/epi2me-labs/pychopper) handles strand orientation and adapter-based read rescue, but direct-cDNA library preparation introduces additional artifact types that Pychopper does not address:

- **Foldback inversions**: the sequenced strand folds back on itself, producing a self-inverted chimeric read.
- **Homopolymer-mediated RT template switching**: during reverse transcription, the RT enzyme detaches at an A/T-rich region on one mRNA and re-primes on another, joining unrelated transcripts into a single chimeric read. These chimeras generate false gene fusion candidates and corrupt isoform quantification.

DirectClean integrates [Breakinator](https://github.com/jheinz27/breakinator) and [Restrander](https://github.com/mritchielab/restrander) with novel detection and rescue algorithms into a single end-to-end pipeline.

### Feature comparison

| Capability | Pychopper | DirectClean |
| :--- | :---: | :---: |
| Strand orientation | ✅ | ✅ |
| Adapter concatemer rescue | ✅ (requires terminal primers) | ✅ (partial internal signal sufficient) |
| Foldback inversion removal | ❌ | ✅ |
| Homopolymer RT template switching detection | ❌ | ✅ |
| Rescue from unclassified reads | ❌ | ✅ |

## Pipeline architecture

| Stage | Name | What it does |
| :--- | :--- | :--- |
| 1 | Breakinator | Remove foldback inversion artifacts |
| 2 | Restrander | Orient reads 5'→3', remove RTP-RTP / TSO-TSO artifacts, set aside unorientable reads |
| 3 | Unknowns Rescue | Recover orientable reads from Restrander unknowns via internal adapter detection and self-orientation |
| 4 | Adapter Rescue | Detect internal TSO/RTP adapters in oriented reads, chop and rescue sub-reads |
| 5 | Homopolymer Rescue | Detect RT template switching at A/T-rich chimeric junctions, chop and rescue sub-reads |

Stages 1–2 remove definitively artifactual or unorientable reads. Stages 3, 4, and 5 never discard reads — they chop chimeric reads at artifact junctions and keep the flanking sub-reads as independent sequences.

### How the homopolymer detector works

After minimap2 splice-aware alignment, DirectClean identifies chimeric reads (those with supplementary alignments mapping to different genomic loci). For each chimeric junction, a 10 bp sliding window scans the flanking sequence on both sides. A junction is flagged as an RT template switching artifact if any window satisfies both criteria:

- A/T base density ≥ 85%
- Longest consecutive A or T run ≥ 5 bp

Flagged reads are chopped at the artifact junction. Sub-reads ≥ 100 bp are written to the output; shorter fragments are discarded. Junctions on non-standard contigs (alt loci, unplaced scaffolds) are excluded via a standard-chromosome whitelist.

## Installation

```bash
# Create environment with all dependencies
mamba env create -f environment.yml
mamba activate directclean

# Install DirectClean
poetry install
```

External tools (minimap2, samtools, breakinator, restrander) are included in the conda environment. To install them separately:

```bash
mamba install -c bioconda minimap2 samtools breakinator
mamba install -c genomedk restrander
```

## Usage

```bash
directclean \
  -i raw_reads.fastq \
  -r genome.fa \
  -o results/ \
  -t 8 \
  -j gencode.v41.bed12
```

The `-j` flag provides a junction BED file for guided alignment (recommended: GENCODE annotation in BED12 format).

### Key parameters

| Flag | Default | Description |
| :--- | :--- | :--- |
| `-i`, `--input` | required | Raw FASTQ from ONT direct-cDNA sequencing |
| `-r`, `--reference` | required | Reference genome FASTA |
| `-o`, `--output` | required | Output directory |
| `-t`, `--threads` | 4 | Threads for minimap2, samtools, breakinator |
| `-j`, `--junc-bed` | none | Junction BED12 for guided alignment |
| `--density-threshold` | 0.85 | A/T density threshold for homopolymer detection |
| `--min-run` | 5 | Minimum consecutive A/T run length |
| `--min-confidence` | 2 | Minimum adapter signals (1–3) required to chop |
| `--context-window` | 50 | Bases flanking each junction for scanning |
| `--html-report` | off | Generate an interactive HTML summary report |

Run `directclean -h` for the full list.

## Output

```text
results/
├── directclean.cleaned.fastq          All clean reads + rescued sub-reads
├── directclean.rescued.fastq          Sub-reads rescued by homopolymer chopping
├── directclean.homopolymer_report.tsv Per-read artifact classification
├── directclean.report.html            Interactive HTML report (if --html-report)
├── intermediates/
│   ├── directclean.no_foldback.fastq       After Stage 1
│   ├── directclean.restranded.fastq        After Stage 2
│   ├── directclean.unknowns_rescued.fastq  Stage 3 output
│   ├── directclean.rescued.fastq           After Stage 4
│   ├── directclean.merged.fastq            Stage 3 + Stage 4 merged
│   └── directclean.aligned.sorted.bam      Minimap2 alignment
└── reports/
    └── directclean.rescue_report.tsv  Stage 4 adapter rescue details
```

The primary output is `directclean.cleaned.fastq`. This file contains all reads that passed the pipeline plus rescued sub-reads from Stages 3, 4, and 5, ready for downstream transcript quantification (e.g., IsoQuant, FLAIR) and gene fusion calling (e.g., FusionSeeker, JAFFAL).

## Citation

If you use DirectClean in your research, please cite our manuscript along with the foundational tools integrated into this pipeline:

- **DirectClean:** Guo, Q., Li, Y., & Yang, R. (2026). DirectClean: a comprehensive preprocessing toolkit for Oxford Nanopore direct-cDNA sequencing. *Manuscript in preparation*.
- **Breakinator:** Heinz, J. M., Meyerson, M., & Li, H. (2026). Detecting foldback artifacts in long-reads. *BMC Genomics*.
- **Restrander:** Schuster, J., Ritchie, M. E., & Gouil, Q. (2023). Restrander: rapid orientation and artefact removal for long-read cDNA data. *NAR Genomics and Bioinformatics*, 5(4), lqad108.

## License

MIT

## Contact

- Qingxiang Guo — qingxiang.guo@northwestern.edu
- Rendong Yang Lab — https://github.com/ylab-hi