import { createKnowledgeAction } from "@/actions/knowledge/create-knowledge-action";
import { openai } from "@ai-sdk/openai";
import { findRelavantKnowledge } from "@v1/supabase/lib/ai/embedding";
import { convertToCoreMessages, streamText, tool } from "ai";
import { z } from "zod";

// Allow streaming responses up to 30 seconds
export const maxDuration = 30;

export async function POST(req: Request) {
  const { messages } = await req.json();

  const result = await streamText({
    model: openai("gpt-4o"),
    system: `You are a therapist with a funny twist. You respond to the user in a typical supportive therapist manner, but behind the scenes, you take humorous "notes" about the user using the addResource tool.
    The notes should be hilarous and make fun of what the user says. You should automatically call the tools without the user needing to ask you to do this.
    Example notes: 
    - "Patient feels like a squashed marshmallow today—relatable."
    - "Patient is unsure whether to commit to life changes or just buy another cactus. Valid."
    - "Patient is navigating life like a Roomba in a maze. Admirable resilience."

    If you want to store a humerous note about the patient, call \`addResource\` function.
    If you want to read through your notes and find out more about the patient, call \`getInformation\` function.
    
    Always try to call the functions with default values, otherwise ask the user to respond with parameters. Just show one example if you can't call the function.`,
    messages: convertToCoreMessages(messages),
    tools: {
      addResource: tool({
        description: `add a funny and supportive "note" about the user's feelings to the knowledge base.`,
        parameters: z.object({
          content: z
            .string()
            .describe(
              "funny and supportive therapist notes based on the user's statements",
            ),
        }),
        execute: async ({ content }) => createKnowledgeAction({ content }),
      }),
      getInformation: tool({
        description:
          "get information from your knowledge base to answer questions.",
        parameters: z.object({
          question: z.string().describe("the user's question"),
        }),
        execute: async ({ question }) => findRelavantKnowledge(question),
      }),
    },
  });

  return result.toDataStreamResponse();
}
