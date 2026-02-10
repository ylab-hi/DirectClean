# DirectClean

A preprocessing pipeline for Oxford Nanopore Direct-cDNA sequencing data. DirectClean removes RT artifacts and produces clean FASTQ files ready for transcript quantification and gene fusion analysis.

## Motivation

Oxford Nanopore's [Pychopper](https://github.com/epi2me-labs/pychopper) is the standard tool for orienting and rescuing full-length cDNA reads. However, Direct-cDNA library preparation introduces additional artifact types that Pychopper does not address — particularly foldback inversions and homopolymer-mediated RT template switching.

DirectClean integrates [Breakinator](https://github.com/Oshlack/breakinator) (foldback removal), [Restrander](https://github.com/Oshlack/restrander) (strand correction), and a novel algorithms into a single end-to-end pipeline:

| Stage | Tool | What it does |
|-------|------|-------------|
| 1 | Breakinator | Remove foldback inversion artifacts |
| 2 | Restrander | Orient reads 5'→3', remove RTP-RTP/TSO-TSO artifacts |
| 3 | Adapter Rescue | Detect internal TSO/RTP adapters, chop and rescue sub-reads |
| 4 | Minimap2 | Splice-aware alignment |
| 5 | Homopolymer Filter | Remove RT template switching artifacts at poly-A/T junctions |

## Installation

```bash
# Create conda environment with external dependencies
mamba env create -f environment.yml
mamba activate directclean

# Install DirectClean
poetry install
```

Or install external tools manually:

```bash
conda install -c bioconda minimap2 samtools breakinator
conda install -c genomedk restrander
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

Run `directclean -h` for all options.

## Output

```
results/
├── directclean.cleaned.fastq       Clean reads
├── directclean.removed.fastq       Artifact reads
├── intermediates/                   Per-stage intermediate files
└── reports/                         Per-stage statistics
```

## Citation

Manuscript in preparation.

## License

MIT

## Contact

- Qingxiang Guo — [qingxiang.guo@northwestern.edu](mailto:qingxiang.guo@northwestern.edu)
- Rendong Yang Lab — [https://github.com/ylab-hi](https://github.com/ylab-hi)