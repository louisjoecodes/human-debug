"use client";

import { Button } from "@v1/ui/button";
import { Input } from "@v1/ui/input";
import { useState } from "react";

export default function VoiceTestPage() {
  const [isPlaying, setIsPlaying] = useState(false);
  const [inputText, setInputText] = useState(
    "Hello! This is a test of the voice streaming API.",
  );

  const handlePlayVoice = async () => {
    setIsPlaying(true);
    try {
      const response = await fetch("/api/voice", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          message: inputText,
          voice: "L0Dsvb3SLTyegXwtm47J",
        }),
      });

      if (!response.ok) {
        throw new Error("Failed to fetch audio stream");
      }

      const audioContext = new window.AudioContext();
      const source = audioContext.createBufferSource();

      const arrayBuffer = await response.arrayBuffer();
      const audioBuffer = await audioContext.decodeAudioData(arrayBuffer);

      source.buffer = audioBuffer;
      source.connect(audioContext.destination);
      source.start(0);

      source.onended = () => {
        setIsPlaying(false);
      };
    } catch (error) {
      console.error("Error playing audio:", error);
      setIsPlaying(false);
    }
  };

  return (
    <div className="flex flex-col items-center justify-center min-h-screen">
      <h1 className="text-2xl font-bold mb-4">Voice Streaming Test</h1>
      <div className="w-full max-w-md mb-4">
        <Input
          value={inputText}
          onChange={(e) => setInputText(e.target.value)}
          placeholder="Enter text to convert to speech"
          className="w-full"
        />
      </div>
      <Button
        onClick={handlePlayVoice}
        disabled={isPlaying || !inputText.trim()}
        className="px-4 py-2 bg-blue-500 text-white rounded"
      >
        {isPlaying ? "Playing..." : "Play Voice"}
      </Button>
    </div>
  );
}
