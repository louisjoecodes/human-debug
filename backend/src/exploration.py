import pysam
import random
from Bio.Seq import Seq

def fetch_sequence_info(chromosome, start, stop):
    with pysam.FastaFile("data/hg38.fa") as reference_genome:
        return reference_genome.fetch(chromosome, start, stop)

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

def generate_cancerous_genome(original_sequence, mutation_type, mutation_position=185):
    """
    Generate a simulated cancerous genome with a specified mutation.
    
    :param original_sequence: The reference BRCA1 sequence (full genome)
    :param mutation_type: The type of mutation to inject (e.g., '185delAG', 'SNP', 'insertion')
    :param mutation_position: The position of the mutation (default is 185)
    :return: The full mutated sequence and the annotation
    """
    # Ensure the mutation position is within the sequence
    if mutation_position >= len(original_sequence):
        raise ValueError("Mutation position is outside the sequence range")
    
    mutated_sequence = list(original_sequence)  # Convert to list for easy manipulation
    annotation = ""

    if mutation_type == '185delAG':
        # Create the mutated sequence with 185delAG mutation
        del mutated_sequence[mutation_position:mutation_position+2]
        annotation = f"{mutation_position}delAG Mutation"
    
    elif mutation_type == 'SNP':
        # Single Nucleotide Polymorphism
        bases = ['A', 'T', 'C', 'G']
        original_base = mutated_sequence[mutation_position]
        new_base = original_base
        while new_base == original_base:
            new_base = random.choice(bases)
        
        mutated_sequence[mutation_position] = new_base
        annotation = f"SNP at position {mutation_position}: {original_base}>{new_base}"
    
    elif mutation_type == 'insertion':
        # Small insertion (e.g., 1-3 bases)
        insertion = ''.join(random.choices(['A', 'T', 'C', 'G'], k=random.randint(1, 3)))
        mutated_sequence[mutation_position:mutation_position] = list(insertion)
        annotation = f"Insertion at position {mutation_position}: {insertion}"
    
    else:
        raise ValueError("Unsupported mutation type")

    return ''.join(mutated_sequence), annotation

if __name__ == "__main__":
    # sub_sequence = fetch_sequence_info(chromosome="chr17", start=43044295, stop=43170245)
    # sub_sequence = fetch_sequence_info(chromosome="chr17", start=43044295, stop=43170245)


    reference_genome = Seq(pysam.FastaFile("data/hg38.fa").fetch("chr17", 43044295, 43170245))
    reference_genome2 = Seq(pysam.FastaFile("data/hg38.fa").fetch("chr17", 43044295, 43170245))

    compare_and_annotate(reference_genome, reference_genome2)
    
    # mutations = [
    #     ('185delAG', 185),
    #     # ('SNP', 1000),
    #     # ('insertion', 5000)
    # ]
    
    # for mutation_type, position in mutations:
    #     cancerous_sequence, annotation = generate_cancerous_genome(sub_sequence, mutation_type, position)
        
    #     print(f"\nMutation: {mutation_type}")
    #     print(f"Length of cancerous sequence: {len(cancerous_sequence)}")
    #     print("Annotation:", annotation)
        
    #     # Compare and annotate
    #     variants = compare_and_annotate(sub_sequence, cancerous_sequence)
        
    #     print("\nDetected Variants:")
    #     for variant in variants:
    #         print(f"Position: {variant['position']}")
    #         print(f"Reference: {variant['reference']}")
    #         print(f"Patient: {variant['patient']}")
    #         print(f"Annotation: {variant['annotation']}")
    #         print("---")

    # # Analyze VCF file
    analyze_vcf_file()

    reference_genome.close()
    reference_genome2.close()