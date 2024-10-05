import os
import asyncio
import json

from dotenv import load_dotenv
import uvicorn
from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware

# Load environment variables from .env file
load_dotenv()

# Initialize the FastAPI app
app = FastAPI(
    title="Chat API with TTS",
    description="An API that streams GPT-4 responses and converts them to speech using ElevenLabs TTS.",
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

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)

def start():
    uvicorn.run("src.main:app", host="0.0.0.0", port=8000, reload=True)

if __name__ == "__main__":
    start()
