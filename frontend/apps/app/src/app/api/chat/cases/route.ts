import { createCaseAction } from "@/actions/cases/create-case-action";
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
    system: `
    You are an administrative assistant in a hospital. You are responsible for creating new cases and getting information about existing cases.

    If you want to create a new case, call \`createCase\` function.
    If you want to get a case, call \`getCase\` function.
    
    Always try to call the functions with default values, otherwise ask the user to respond with parameters. Just show one example if you can't call the function.`,
    messages: convertToCoreMessages(messages),
    tools: {
      createCase: tool({
        description: `create a new case in the database.`,
        parameters: z.object({
          first_name: z.string().describe("The first name of the patient"),
          last_name: z.string().describe("The last name of the patient"),
          date_of_birth: z
            .string()
            .describe("The date of birth of the patient"),
        }),
        execute: async ({ first_name, last_name, date_of_birth }) =>
          createCaseAction({ first_name, last_name, date_of_birth }),
      }),
      // getInformation: tool({
      //   description:
      //     "get information from your knowledge base to answer questions.",
      //   parameters: z.object({
      //     question: z.string().describe("the user's question"),
      //   }),
      //   execute: async ({ question }) => findRelavantKnowledge(question),
      // }),
    },
  });

  return result.toDataStreamResponse();
}
