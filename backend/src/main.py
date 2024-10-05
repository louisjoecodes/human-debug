import os

import base64
from fastapi import UploadFile, File
from mistralai import Mistral


import utils

from dotenv import load_dotenv
import uvicorn
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

@app.get("/extract_letter_content")
async def get_extract_letter_content():
    return {"message": "Please use POST to upload a file"}

@app.post("/extract_letter_content")
async def extract_letter_content(file: UploadFile = File(...)):
    # content = await file.read()

    image_path = "data/doctor_letter_scan_0.jpg"
    print(f"image_path: {image_path}")

    # Getting the base64 string
    encoded_content = utils.encode_image(image_path)
    print(f"encoded_content: {encoded_content}")
    
    # encoded_content = base64.b64encode(content).decode('utf-8')
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

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)

def start():
    uvicorn.run("src.main:app", host="0.0.0.0", port=8000, reload=True)

if __name__ == "__main__":
    start()
