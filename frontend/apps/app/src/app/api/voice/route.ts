import { env } from "@/env.mjs";
import { ElevenLabsClient, ElevenLabsError } from "elevenlabs";

export async function POST(req: Request) {
  const { message, voice } = await req.json();

  const elevenlabs = new ElevenLabsClient({
    apiKey: env.ELEVENLABS_API_KEY,
  });

  try {
    const audioStream = await elevenlabs.generate({
      voice,
      model_id: "eleven_turbo_v2",
      voice_settings: { similarity_boost: 0.5, stability: 0.5 },
      text: message,
      stream: true, // Enable streaming
    });

    // Create a ReadableStream from the audio stream
    const stream = new ReadableStream({
      async start(controller) {
        for await (const chunk of audioStream) {
          controller.enqueue(chunk);
        }
        controller.close();
      },
    });

    // Return the streaming response
    return new Response(stream, {
      headers: {
        "Content-Type": "audio/mpeg",
        "Transfer-Encoding": "chunked",
      },
    });
  } catch (error) {
    if (error instanceof ElevenLabsError) {
      console.error(error);
      return Response.json(error, { status: error.statusCode });
    }
    console.error(error);
    return Response.json(error, { status: 500 });
  }
}
