import os

from fastapi import UploadFile, File
from mistralai import Mistral
import instructor
import base64
from typing import Literal

from dotenv import load_dotenv
import uvicorn
from pydantic import BaseModel, Field
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

load_dotenv()
client = Mistral(api_key=os.getenv("MISTRAL_API_KEY"))

app = FastAPI(
    title="AI-Assisted Genetic Diagnostics API",
    description="An API that leverages AI to assist doctors in genetic diagnostics and interpretation of medical data.",
    version="1.0.0"
)

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
async def process_letter():

    content = extract_letter_content()
    patient = extract_patient_info(content)

    print(f"content: {content}")
    print(f"patient: {patient.model_dump()}")


    return patient



@app.post("/structure_medical_history")
async def structure_medical_history(content: str):

    return {"message": "Please use POST to send in text"}

async def extract_letter_content(file=None):

    image_path = "^data/doctor_letter_scan_0.jpg"
    encoded_content = encode_image(image_path)

    messages = [
        {
            "role": "user",
            "content": [
                {
                    "type": "text",
                    "text": "Please transcribe the content of this image or PDF."
                },
                {
                    "type": "image_url",
                    "image_url": f"data:image/jpeg;base64,{encoded_content}"
                }
            ]
        }
    ]
    
    chat_response = client.chat.complete(
        model="pixtral-12b-2409",
        messages=messages
    )
    
    transcribed_content = chat_response.choices[0].message.content
    print(transcribed_content)
    
    return {"content": transcribed_content}

def extract_patient_info(content):

    client = instructor.from_mistral(Mistral(api_key=os.getenv("MISTRAL_API_KEY")))

    with open("data/doctor_letter.txt", "r") as file:
        content = file.read()

    class Patient(BaseModel):
        first_name: str
        last_name: str
        date_of_birth: str
        gender: Literal["Female", "Male"]
        age: int
        disease: str = Field(..., description="The main reason for the genomic test")
        

    patient = client.chat.completions.create(
        model="pixtral-12b-2409",
        response_model=Patient,
        messages=[
            {"role": "user", "content": f"Extract patient information from this text:\n\n{content}"}
        ],
    )

    print(patient.model_dump())

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

def start():
    uvicorn.run("src.main:app", host="0.0.0.0", port=8000, reload=True)

if __name__ == "__main__":
    start()
