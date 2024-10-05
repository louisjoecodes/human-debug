#!/bin/bash

# Input FASTQ file
FASTQ_FILE="SRR765993_2.fastq.gz"

# Reference genome (you need to provide the path to your reference genome)
REFERENCE_GENOME="data/hg38.fa"

# Output BAM file name
OUTPUT_BAM="${FASTQ_FILE%.fastq.gz}.bam"

# Step 1: Align reads to the reference genome using BWA-MEM
bwa mem $REFERENCE_GENOME $FASTQ_FILE > aligned_reads.sam

# Step 2: Convert SAM to BAM
samtools view -bS aligned_reads.sam > aligned_reads.bam

# Step 3: Sort the BAM file
samtools sort aligned_reads.bam -o sorted_aligned_reads.bam

# Step 4: Index the sorted BAM file
samtools index sorted_aligned_reads.bam

# The final BAM file is now ready
mv sorted_aligned_reads.bam $OUTPUT_BAM
mv sorted_aligned_reads.bam.bai ${OUTPUT_BAM}.bai

# Clean up intermediate files
rm aligned_reads.sam aligned_reads.bam

echo "Conversion complete. Output BAM file: $OUTPUT_BAM"