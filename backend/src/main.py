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
from fastapi import FastAPI, File, UploadFile
from fastapi.middleware.cors import CORSMiddleware

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
    # comment: str | None = Field(None, description="Additional comments about the phenotype")
    # descendant_count: int = Field(..., description="The number of descendants for this phenotype")
    # synonyms: List[str] = Field(default_factory=list, description="Alternative names for the phenotype")
    # xrefs: List[str] = Field(default_factory=list, description="Cross-references to other databases")
    # translations: Any | None = Field(None, description="Translations of the phenotype name, if available")

class Disease(BaseModel):
    id: str = Field(..., description="The identifier for the disease")
    name: str = Field(..., description="The name of the disease")
    mondoId: str = Field(..., description="The MONDO identifier for the disease")
    description: str | None = Field(None, description="Optional description of the disease")


class Gene(BaseModel):
    id: str = Field(..., description="The identifier for the gene, e.g., 'NCBIGene:4340'")
    name: str = Field(..., description="The name of the gene, e.g., 'MOG'")


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
    return {"message": "Hello World"}

@app.get("/process_letter")
async def get_process_letter():
    return {"message": "Please use POST to upload a file"}

@app.post("/process_letter")
async def process_letter(file: UploadFile):
    print("PROCESSING LETTER", file)
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

    with open("data/doctor_letter.txt", "r") as file:
        content = file.read()

    patient = client.chat.completions.create(
        model="pixtral-12b-2409",
        response_model=Patient,
        messages=[
            {"role": "user", "content": f"Extract patient information from this text:\n\n{content}"}
        ],
    )

    return patient


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

def start():
    uvicorn.run("src.main:app", host="0.0.0.0", port=8000, reload=True)

if __name__ == "__main__":
    start()
