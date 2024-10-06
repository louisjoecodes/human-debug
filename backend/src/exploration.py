import pysam
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess
import sys

def fetch_genomic_reference_sequence(chromosome, start=None, stop=None) -> Seq:
    with pysam.FastaFile("data/hg38.fa") as reference_genome:
        start = start or 0
        stop = stop or reference_genome.get_reference_length(chromosome)
        return Seq(reference_genome.fetch(chromosome, start, stop))

def create_fasta(sequence, filename, seq_id):
    record = SeqRecord(sequence, id=seq_id, description="")
    SeqIO.write(record, filename, "fasta")

    
def create_mutated_sequence(sequence, mutation_position, deletion_length):
    return sequence[:mutation_position] + sequence[mutation_position + deletion_length:]

def align_sequence(reference_file, query_file, bam_file):
    """Align sequence from query file to reference and create a BAM file using BWA."""
    
    try:
        # Index the reference file (if not already indexed)
        subprocess.run(["bwa", "index", reference_file], check=True)
        
        # Perform alignment using BWA
        sam_file = bam_file.replace('.bam', '.sam')
        with open(sam_file, 'w') as sam_output:
            subprocess.run(["bwa", "mem", reference_file, query_file], stdout=sam_output, check=True)
        
        # Convert SAM to BAM
        pysam.view("-b", "-o", bam_file, sam_file, catch_stdout=False)
        
        # Index the BAM file
        pysam.index(bam_file)
        
        # Clean up the intermediate SAM file
        subprocess.run(["rm", sam_file], check=True)

        print(f"Alignment completed. BAM file created: {bam_file}")
    
    except FileNotFoundError:
        print("Error: BWA not found. Please install BWA and ensure it's in your system PATH.")
        print("Installation instructions:")
        print("  - macOS (using Homebrew): brew install bwa")
        print("  - Ubuntu/Debian: sudo apt-get install bwa")
        print("  - From source: Visit http://bio-bwa.sourceforge.net/ for instructions")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error during alignment process: {e}")
        sys.exit(1)

def analyze_vcf_file():
    with pysam.VariantFile("data/clinvar.vcf.gz") as vcf:
        print(f"Number of samples: {len(vcf.header.samples)}")
        
        # Print the first few variants
        print("\nFirst 5 variants:")
        for i, record in enumerate(vcf.fetch()):
            if i >= 5:
                break
            print(f"Chromosome: {record.chrom}, Position: {record.pos}, ID: {record.id}")
            print(f"Reference: {record.ref}, Alternative: {', '.join(record.alts)}")
            print(f"CLNSIG: {record.info.get('CLNSIG', 'Not available')}")
            print("---")

        # Example: Count variants by chromosome
        chrom_counts = {}
        vcf.reset()
        for record in vcf.fetch():
            chrom_counts[record.chrom] = chrom_counts.get(record.chrom, 0) + 1

        print("\nVariant counts by chromosome:")
        for chrom, count in chrom_counts.items():
            print(f"{chrom}: {count}")

    return chrom_counts

def index_reference(reference_file):
    """Index the reference genome if it's not already indexed."""
    try:
        subprocess.run(["samtools", "faidx", reference_file], check=True)
    except subprocess.CalledProcessError:
        print(f"Error indexing reference file: {reference_file}")
        raise

def call_variants(reference_file, bam_file, output_vcf):
    """Perform variant calling using FreeBayes."""
    try:
        # Ensure the reference is indexed
        index_reference(reference_file)
        
        # Run FreeBayes with less stringent parameters
        command = [
            "freebayes",
            "-f", reference_file,
            "-C", "1",  # Minimum alternate count (reduced from 5)
            "--min-alternate-fraction", "0.01",  # Reduced from 0.05
            "--min-base-quality", "10",
            "--min-mapping-quality", "10",
            "-b", bam_file,
            "-v", output_vcf
        ]
        
        subprocess.run(command, check=True)
        print(f"Variant calling completed. VCF file created: {output_vcf}")
        
    except subprocess.CalledProcessError as e:
        print(f"Error during variant calling: {e}")
        raise

def filter_variants(input_vcf, output_filtered_vcf):
    """Filter variants based on quality and depth."""
    try:
        command = [
            "bcftools", "filter",
            "-i", "QUAL>20 && DP>10",
            "-o", output_filtered_vcf,
            "-O", "v",
            input_vcf
        ]
        
        subprocess.run(command, check=True)
        print(f"Variant filtering completed. Filtered VCF file created: {output_filtered_vcf}")
        
    except subprocess.CalledProcessError as e:
        print(f"Error during variant filtering: {e}")
        raise

def analyze_variants(vcf_file):
    """Basic analysis of the variants in the VCF file."""
    try:
        with pysam.VariantFile(vcf_file) as vcf:
            variant_count = 0
            snp_count = 0
            indel_count = 0
            
            for record in vcf:
                variant_count += 1
                if len(record.ref) == 1 and len(record.alts[0]) == 1:
                    snp_count += 1
                else:
                    indel_count += 1
            
            print(f"Total variants: {variant_count}")
            print(f"SNPs: {snp_count}")
            print(f"Indels: {indel_count}")
    
    except Exception as e:
        print(f"Error during variant analysis: {e}")
        raise


def introduce_random_mutations(sequence, max_mutations=5):
    """
    Introduce a limited number of random mutations (noise) into the given sequence.
    
    :param sequence: The input DNA sequence (Seq object)
    :param max_mutations: The maximum number of mutations to introduce
    :return: A new Seq object with random mutations
    """
    nucleotides = ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']
    mutated_sequence = list(str(sequence))
    sequence_length = len(mutated_sequence)
    
    for _ in range(max_mutations):
        position = random.randint(0, sequence_length - 1)
        mutation_type = random.choice(['add', 'delete', 'switch'])
        
        if mutation_type == 'add':
            mutated_sequence.insert(position, random.choice(nucleotides))
            sequence_length += 1
        elif mutation_type == 'delete':
            del mutated_sequence[position]
            sequence_length -= 1
        else:  # switch
            mutated_sequence[position] = random.choice(nucleotides)
    
    return Seq(''.join(mutated_sequence))


if __name__ == "__main__":

    reference_file = "data/reference_brca1.fa"
    patient_file = "data/patient_brca1.fa"
    bam_file = "data/aligned_brca1.bam"
    output_vcf = "data/variants.vcf"
    filtered_vcf = "data/filtered_variants.vcf"

    reference_sequence_brca1 = fetch_genomic_reference_sequence(chromosome="chr17", start=200000, stop=600000) #start=43044295, stop=43170245
    create_fasta(reference_sequence_brca1, reference_file, "reference")

    # print(reference_sequence_brca1)

    # mutated_sequence_brca1 = create_mutated_sequence(reference_sequence_brca1, mutation_position=185, deletion_length=2)
    noisy_mutated_sequence_brca1 = introduce_random_mutations(reference_sequence_brca1, max_mutations=1)

    print("Sequences are identical:", noisy_mutated_sequence_brca1 == reference_sequence_brca1)
    if noisy_mutated_sequence_brca1 != reference_sequence_brca1:
        print("Number of differences:", sum(1 for a, b in zip(str(noisy_mutated_sequence_brca1), str(reference_sequence_brca1)) if a != b))
    create_fasta(noisy_mutated_sequence_brca1, patient_file, "patient")

    # print("--------------------------")
    # print(noisy_mutated_sequence_brca1)
    align_sequence(reference_file=reference_file, query_file=patient_file, bam_file=bam_file)


    def check_bam_file(bam_file):
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            print(f"Checking BAM file: {bam_file}")
            print(f"Number of mapped reads: {bam.mapped}")
            print(f"Number of unmapped reads: {bam.unmapped}")
            
            mismatches = 0
            for read in bam.fetch():
                if read.has_tag("NM"):
                    mismatches += read.get_tag("NM")
        print(f"Total mismatches: {mismatches}")

    # Add this to your main function
    check_bam_file(bam_file)

    # New variant calling, filtering, and analysis workflow
    call_variants(reference_file, bam_file, output_vcf)
    # filter_variants(output_vcf, filtered_vcf)
    # analyze_variants(filtered_vcf)
