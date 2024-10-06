import os
import requests
from mistralai import Mistral
import instructor
import base64
import pprint
import io
from pdf2image import convert_from_bytes
from typing import Literal
from urllib.parse import quote

from PIL import Image
from dotenv import load_dotenv
import uvicorn
from typing import Dict, List, Any
from pydantic import BaseModel, Field
from fastapi import FastAPI, UploadFile
from fastapi.middleware.cors import CORSMiddleware

import subprocess
import vcfpy

load_dotenv()
pp = pprint.PrettyPrinter(indent=4)

client = Mistral(api_key=os.getenv("MISTRAL_API_KEY"))

app = FastAPI(
    title="AI-Assisted Genetic Diagnostics API",
    description="An API that leverages AI to assist doctors in genetic diagnostics and interpretation of medical data.",
    version="1.0.0"
)


class Patient(BaseModel):
    first_name: str
    last_name: str
    date_of_birth: str = Field(..., pattern=r'^\d{4}-\d{2}-\d{2}$', description="Date of birth in YYYY-MM-DD format")
    gender: Literal["Female", "Male"]
    age: int
    disease: str = Field(..., description="The disease for ordering genomic test. Max 1-2 words")

class Phenotype(BaseModel):
    id: str = Field(..., description="The identifier for the phenotype")
    name: str = Field(..., description="The name of the phenotype")
    definition: str | None = Field(None, description="The definition of the phenotype")
    comment: str | None = Field(None, description="Additional comments about the phenotype")
    descendant_count: int = Field(..., description="The number of descendants for this phenotype")
    synonyms: List[str] = Field(default_factory=list, description="Alternative names for the phenotype")
    xrefs: List[str] = Field(default_factory=list, description="Cross-references to other databases")
    translations: Any | None = Field(None, description="Translations of the phenotype name, if available")

class Disease(BaseModel):
    id: str = Field(..., description="The identifier for the disease")
    name: str = Field(..., description="The name of the disease")
    mondoId: str = Field(..., description="The MONDO identifier for the disease")
    description: str | None = Field(None, description="Optional description of the disease")


class Gene(BaseModel):
    id: str = Field(..., description="The identifier for the gene, e.g., 'NCBIGene:4340'")
    name: str = Field(..., description="The name of the gene, e.g., 'MOG'")


class Variant(BaseModel):
    chromosome: str = Field(..., description="Chromosome where the variant is located")
    position: int = Field(..., description="Position of the variant on the chromosome")
    reference: str = Field(..., description="Reference allele")
    alternate: str = Field(..., description="Alternate allele")
    gene: str = Field(..., description="Gene affected by the variant")
    consequence: str = Field(..., description="Predicted consequence of the variant")
    significance: Literal[
        "pathogenic",
        "likely_pathogenic",
        "uncertain_significance",
        "likely_benign",
        "benign",
        "not_provided",
        "drug_response",
        "conflicting",
        "risk_factor",
        "association",
        "protective"
    ] = Field(..., description="Clinical significance of the variant")


# Allow all CORS origins (for demo purposes only)
app.add_middleware(
    CORSMiddleware,
    allow_credentials=True,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
async def root():
    return {"status": "ok"}

@app.get("/process_letter")
async def get_process_letter():
    return {"message": "Please use POST to upload a file"}

@app.post("/process_letter")
async def process_letter(file: UploadFile):
    content = extract_letter_content(file)
    patient = extract_patient_info(content)
    phenotype_classes = fetch_phenotype_classes(patient.disease)

    print(f"content: {content}")
    print(f"patient: {patient.model_dump()}")
    print(f"phenotype_classes: {phenotype_classes}")
    
    return {
        
        "patient": patient.model_dump(),
        "phenotype_classes": phenotype_classes
    }

@app.post("/extract_text")
async def extract_text(file: UploadFile):
    content = await extract_letter_content(file)
    return content

@app.post("/analyze")
async def analyze():
    phenotype_ids = ["HP:0003002", "HP:0010619"]
    annotations = list(map(lambda phenotype_id: fetch_annotation(phenotype_id=phenotype_id), phenotype_ids))

    pp.pprint(annotations)
    return annotations

@app.get("/gene_to_phenotypes/{gene_id}")
async def gene_to_phenotypes(gene_id: str="NCBIGene:3161"):
    """
    Fetch associated phenotypes and diseases for a given gene ID.
    
    Args:
    gene_id (str): The NCBI Gene ID, e.g., "NCBIGene:3161"

    Returns:
    Dict[str, List[Dict[str, Any]]]: A dictionary containing lists of associated diseases and phenotypes
    """
    base_url = "https://ontology.jax.org/api/network/annotation"
    encoded_gene_id = quote(gene_id)
    url = f"{base_url}/{encoded_gene_id}"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        
        return {
            "diseases": [Disease(**disease) for disease in data.get("diseases", [])],
            "phenotypes": [Phenotype(**phenotype) for phenotype in data.get("phenotypes", [])]
        }
    except requests.RequestException as e:
        print(f"Error fetching gene: {e}")
        return {"diseases": [], "phenotypes": []}



async def extract_letter_content(file: UploadFile):
    content = await file.read()
    images: List[Image.Image] = []

    if file.content_type == "application/pdf":
        images = convert_from_bytes(content, fmt='jpeg')
    elif file.content_type in ["image/jpeg", "image/png"]:
        images = [Image.open(io.BytesIO(content))]
    else:
        raise ValueError(f"Unsupported file type: {file.content_type}")

    transcribed_contents: List[str] = []

    for img in images:
        img_byte_arr = io.BytesIO()
        img.save(img_byte_arr, format='JPEG')
        encoded_content = base64.b64encode(img_byte_arr.getvalue()).decode('utf-8')

        message = {
            "role": "user",
            "content": [
                {"type": "text", "text": "Please transcribe the content of this image."},
                {"type": "image_url", "image_url": f"data:image/jpeg;base64,{encoded_content}"}
            ]
        }
        
        chat_response = client.chat.complete(
            model="pixtral-12b-2409",
            messages=[message]
        )
        
        transcribed_contents.append(chat_response.choices[0].message.content)
    
    full_transcription = "\n\n".join(transcribed_contents)
    return {"content": full_transcription}

def extract_patient_info(content):
    client = instructor.from_mistral(Mistral(api_key=os.getenv("MISTRAL_API_KEY")))
    return client.chat.completions.create(
        model="pixtral-12b-2409",
        response_model=Patient,
        messages=[
            {"role": "user", "content": f"Extract patient information from this text, her name is always 'Sarah Elizabeth Thompson' and her DOB is always '03/15/1978':\n\n{content}"}
        ],
    )


def encode_image(image_path):
    """Encode the image to base64."""
    try:
        with open(image_path, "rb") as image_file:
            return base64.b64encode(image_file.read()).decode('utf-8')
    except FileNotFoundError:
        print(f"Error: The file {image_path} was not found.")
        return None
    except Exception as e:
        print(f"Error: {e}")
        return None

def fetch_phenotype_classes(disease):
    """Fetch related terms from the JAX ontology based on the patient's disease."""
    base_url = "https://ontology.jax.org/api/hp/search"
    encoded_disease = quote(disease)
    url = f"{base_url}?q={encoded_disease}&page=1&limit=10"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        terms = data.get("terms", [])
        
        # Extract relevant information from each term
        simplified_terms = [
            {
                "id": term["id"],
                "name": term["name"],
                "definition": term["definition"]
            }
            for term in terms
        ]
        
        return simplified_terms
    except requests.RequestException as e:
        print(f"Error fetching ontology terms: {e}")
        return []
    

def fetch_annotation(phenotype_id: str="HP:0003002") -> Dict[str, List[Dict[str, Any]]]:
    """
    Fetch related diseases, genes, and medical actions for a given phenotype ID from the JAX ontology.
    
    Args:
    phenotype_id (str): The HP (Human Phenotype) ID to query, e.g., "HP:0003002"

    Returns:
    Dict[str, List[Dict[str, Any]]]: A dictionary containing lists of related diseases, genes, and medical actions
    """
    base_url = "https://ontology.jax.org/api/network/annotation"
    encoded_phenotype_id = quote(phenotype_id)
    url = f"{base_url}/{encoded_phenotype_id}"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        
        return {
            "diseases": data.get("diseases", []),
            "genes": data.get("genes", []),
            "medical_actions": data.get("medicalActions", [])
        }
    except requests.RequestException as e:
        print(f"Error fetching phenotype network: {e}")
        return {"diseases": [], "genes": [], "medical_actions": []}

@app.get("/variants")
async def get_variants() -> List[Variant]:
    """
    Return a list of mocked variants identified from genomic analysis.
    This endpoint is for frontend development purposes.
    """
    mock_variants = [
        Variant(
            chromosome="1",
            position=69897,
            reference="A",
            alternate="G",
            gene="OR4F5",
            consequence="missense_variant",
            significance="uncertain_significance"
        ),
        Variant(
            chromosome="7",
            position=117199644,
            reference="ATCT",
            alternate="A",
            gene="CFTR",
            consequence="frameshift_variant",
            significance="pathogenic"
        ),
        Variant(
            chromosome="17",
            position=41197708,
            reference="G",
            alternate="A",
            gene="BRCA1",
            consequence="stop_gained",
            significance="likely_pathogenic"
        ),
        Variant(
            chromosome="X",
            position=31496081,
            reference="C",
            alternate="T",
            gene="DMD",
            consequence="splice_donor_variant",
            significance="pathogenic"
        ),
        Variant(
            chromosome="4",
            position=1807894,
            reference="G",
            alternate="A",
            gene="FGFR3",
            consequence="missense_variant",
            significance="benign"
        )
    ]
    return mock_variants

@app.post("/analyze_gene_sequence")
async def analyze_gene_sequence(upload_file: UploadFile):
    """Analyze the predefined .fastq gene sequence file and return identified variants."""

    # query_file = "data/Run2_IonXpress_008.fastq"
    reference_file = "data/hg38.fa"
    
    bam_file = align_sequence(query_file=upload_file, reference_file=reference_file)
    print(f"BAM file created: {bam_file}")

    vcf_file = call_variants(input_bam=bam_file, reference_file=reference_file)
    print(f"VCF file created: {vcf_file}")

    filtered_vcf_file = filter_variants(input_vcf=vcf_file)
    print(f"Filtered VCF file created: {filtered_vcf_file}")

    variants = parse_variants(filtered_vcf_file)
    print(f"Variants parsed: {[variant.model_dump() for variant in variants]}")

    return {"variants": [variant.model_dump() for variant in variants]}

def align_sequence(query_file, reference_file="data/hg38.fa"):
    """Align sequence from FASTQ query file to reference and create a BAM file using BWA."""

    return "data/output.sorted.bam"
    
    # TODO: Reactivate this code for inference
    # bam_file = f"data/{os.path.splitext(os.path.basename(query_file))[0]}.bam"
    # try:
    #     # Index the reference file if not already indexed
    #     # if not os.path.exists(reference_file + ".bwt"):
    #     #     print(f"Indexing reference file: {reference_file}")
    #     #     subprocess.run(["bwa", "index", reference_file], check=True)

    #     # if not os.path.exists(reference_file + ".fa.fai"):
    #     #     print(f"Indexing reference file: {reference_file}")
    #     #     subprocess.run(["bwa", "index", reference_file], check=True)  

    #     # Perform alignment using BWA
    #     sam_file = bam_file.replace('.bam', '.sam')
    #     print(f"Aligning sequences to reference: {query_file} -> {reference_file}")
    #     with open(sam_file, 'w') as sam_output:
    #         result = subprocess.run([
    #             "bwa", "mem",
    #             "-A", "1",  # Match score
    #             "-B", "2",  # Mismatch penalty
    #             "-O", "3",  # Gap open penalty
    #             "-E", "1",  # Gap extension penalty
    #             "-L", "0",  # Clipping penalty
    #             "-T", "10", # Minimum score to output
    #             "-a",       # Output all alignments for SE or unpaired PE
    #             reference_file, query_file
    #         ], stdout=sam_output, stderr=subprocess.PIPE, text=True, check=True)
        
    #     print(f"BWA mem stdout: {result.stdout}")
    #     print(f"BWA mem stderr: {result.stderr}")

    #     # Convert SAM to BAM
    #     print(f"Converting SAM to BAM: {sam_file} -> {bam_file}")
    #     pysam.view("-b", "-o", bam_file, sam_file, catch_stdout=False)
        
    #     # Sort the BAM file
    #     sorted_bam = bam_file.replace('.bam', '.sorted.bam')
    #     print(f"Sorting BAM file: {bam_file} -> {sorted_bam}")
    #     pysam.sort("-o", sorted_bam, bam_file)
        
    #     # Replace the original BAM file with the sorted one
    #     subprocess.run(["mv", sorted_bam, bam_file], check=True)
        
    #     # Index the sorted BAM file
    #     print(f"Indexing BAM file: {bam_file}")
    #     pysam.index(bam_file)
        
    #     # Clean up the intermediate SAM file
    #     subprocess.run(["rm", sam_file], check=True)

    #     print(f"Alignment completed. BAM file created: {bam_file}")
        # return bam_file

    # except subprocess.CalledProcessError as e:
    #     print(f"Command failed: {' '.join(e.cmd)}", file=sys.stderr)
    #     print(f"Exit status: {e.returncode}", file=sys.stderr)
    #     print(f"STDOUT: {e.stdout}", file=sys.stderr)
    #     print(f"STDERR: {e.stderr}", file=sys.stderr)
    #     raise
    # except Exception as e:
    #     print(f"An unexpected error occurred: {str(e)}", file=sys.stderr)
    #     raise

def call_variants(input_bam, reference_file, output_vcf=None) -> str:
    """Perform variant calling using FreeBayes on a BAM file with a reference genome."""
    try:
        # Generate output_vcf based on the input BAM file name if not provided
        if output_vcf is None:
            bam_basename = os.path.splitext(os.path.basename(input_bam))[0]
            output_vcf = f"data/variants_{bam_basename}.vcf"

        # Run FreeBayes with less stringent parameters
        command = [
            "freebayes",
            "-f", reference_file,
            "-C", "1",
            "--min-alternate-fraction", "0.005",
            "--min-base-quality", "5",
            "--min-mapping-quality", "5",
            "--min-coverage", "1",
            "--min-alternate-count", "1",
            "--haplotype-length", "0",
            "--pooled-continuous",
            "--use-best-n-alleles", "5",
            "-b", input_bam,
            "-v", output_vcf
        ]
        subprocess.run(command, check=True)
        print(f"Variant calling completed. VCF file created: {output_vcf}")
        return output_vcf
        
    except subprocess.CalledProcessError as e:
        print(f"Error during variant calling: {e}")
        raise

def filter_variants(input_vcf: str) -> str:
    """Filter variants based on quality."""
    try:
        base_name = os.path.splitext(input_vcf)[0]
        output_filtered_vcf = f"{base_name}_filtered.vcf"

        command = [
            "bcftools", "filter",
            "-i", "QUAL>20",
            "-o", output_filtered_vcf,
            "-O", "v",
            input_vcf
        ]
        
        subprocess.run(command, check=True)
        print(f"Variant filtering completed. Filtered VCF file created: {output_filtered_vcf}")
        
        return output_filtered_vcf
        
    except subprocess.CalledProcessError as e:
        print(f"Error during variant filtering: {e}")
        raise

def parse_variants(vcf_file: str) -> List[Variant]:
    """
    Parse a VCF file and return a list of Variant objects.
    """
    variants = []
    reader = vcfpy.Reader.from_path(vcf_file)
    
    for record in reader:
        # Extract basic information
        chromosome = record.CHROM
        position = record.POS
        reference = record.REF
        alternate = ','.join(str(alt.value) for alt in record.ALT)
        
        # Extract variant information
        gene = record.INFO.get('GENE', "Unknown")
        consequence = record.INFO.get('CONSEQUENCE', "Unknown")
        
        # Determine clinical significance
        significance = "uncertain_significance"
        if 'CLNSIG' in record.INFO:
            significance = record.INFO['CLNSIG'].lower()

        variants.append(Variant(
            chromosome=chromosome,
            position=position,
            reference=reference,
            alternate=alternate,
            gene=gene,
            consequence=consequence,
            significance=significance
        ))
    
    return variants

def start():
    uvicorn.run("src.main:app", host="0.0.0.0", port=8000, reload=True)

if __name__ == "__main__":
    start()