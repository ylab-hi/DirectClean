"""Smoke test — run the full DirectClean pipeline on mini data."""
import sys
import logging

sys.path.insert(0, "../src")

from directclean.pipeline import DirectCleanPipeline, PipelineConfig

logging.basicConfig(level=logging.INFO, format="%(message)s")

# === Paths — edit these ===
INPUT_FASTQ = "mini.fastq"
REFERENCE = "/gpfs/projects/b1171/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa"
JUNC_BED = "/gpfs/projects/b1171/qgn1237/1_my_database/GRCh38_p13/gencode.v41.bed12"
OUTPUT_DIR = "test_output"

# === Run full pipeline ===
config = PipelineConfig(
    threads=8,
    junc_bed=JUNC_BED,
)

pipeline = DirectCleanPipeline(
    input_fastq=INPUT_FASTQ,
    reference=REFERENCE,
    output_dir=OUTPUT_DIR,
    config=config,
)
report = pipeline.run()
print(report)