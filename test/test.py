import sys
sys.path.insert(0, "../src")

from directclean.filter import ArtifactClassifier, HomopolymerConfig

classifier = ArtifactClassifier(
    bam_path="mini.bam",
    input_fastq="mini.fastq",
    output_dir="test_output/",
)
report = classifier.run()
print(report)
