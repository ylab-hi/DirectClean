"""Smoke test — run the full DirectClean pipeline on mini data."""

import sys
import logging

# If not installed via poetry, fallback to src path
sys.path.insert(0, "../src")

from directclean.pipeline import DirectCleanPipeline, PipelineConfig
from directclean.rescuer.adaptor_seq import AdapterConfig
from directclean.filter.homopolymer import HomopolymerConfig

logging.basicConfig(level=logging.INFO, format="%(message)s")

# === Paths — edit these ===
INPUT_FASTQ = "mini.fastq"
REFERENCE = "/gpfs/projects/b1171/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa" 
OUTPUT_DIR = "test_output"

# === Run full pipeline ===
pipeline = DirectCleanPipeline(
    input_fastq=INPUT_FASTQ,
    reference=REFERENCE,
    output_dir=OUTPUT_DIR,
    threads=8,
)
report = pipeline.run()
print(report)
