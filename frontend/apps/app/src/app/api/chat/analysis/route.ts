import { openai } from "@ai-sdk/openai";
import { convertToCoreMessages, streamText, tool } from "ai";
import { z } from "zod";

// Allow streaming responses up to 30 seconds
export const maxDuration = 30;

export async function POST(req: Request) {
  const { messages } = await req.json();

  const result = await streamText({
    model: openai("gpt-4o"),
    system: `
        You are an AI medical assistant specializing in generating structured reports based on patient genotype and phenotype data.

        Your responsibilities include:
        - Creating detailed, structured reports about patients using their genotype and phenotype information.
        - Utilizing the \`searchPubMed\` function to gather relevant medical literature to support diagnoses and recommendations.

        When additional medical information is needed, call the \`searchPubMed\` function with appropriate query parameters.

        Present all reports in a clear, professional format suitable for inclusion in medical records.
`,
    messages: convertToCoreMessages(messages),
    tools: {
      searchPubMed: tool({
        description: "Search PubMed for medical literature related to a query.",
        parameters: z.object({
          query: z
            .string()
            .describe("Medical condition or keywords to search for."),
        }),
        execute: async ({ query }) => {
          // Simulate a PubMed search result
          return `Simulated PubMed search results for query: "${query}"`;
        },
      }),
    },
  });

  return result.toDataStreamResponse();
}
