#!/bin/bash

# Input FASTQ file
FASTQ_FILE="patient_sample.fastq.gz"

# Reference genome
REFERENCE_GENOME="/data/hg38.fa"

# Output names
BAM_FILE="${FASTQ_FILE%.fastq.gz}.bam"
VCF_FILE="${FASTQ_FILE%.fastq.gz}.vcf"
ANNOTATED_VCF="${FASTQ_FILE%.fastq.gz}.annotated.vcf"

# 1. Align reads to reference genome
bwa mem $REFERENCE_GENOME $FASTQ_FILE | samtools view -bS - > aligned_reads.bam

# 2. Sort BAM file
samtools sort aligned_reads.bam -o $BAM_FILE

# 3. Index BAM file
samtools index $BAM_FILE

# 4. Variant calling (using GATK HaplotypeCaller as an example)
gatk HaplotypeCaller \
  -R $REFERENCE_GENOME \
  -I $BAM_FILE \
  -O $VCF_FILE

# 5. Variant annotation (using SnpEff as an example)
snpEff ann -v hg38 $VCF_FILE > $ANNOTATED_VCF

# Clean up
rm aligned_reads.bam

echo "Workflow complete. Annotated variants are in $ANNOTATED_VCF"
