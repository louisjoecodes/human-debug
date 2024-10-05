import { openai } from "@ai-sdk/openai";
import { createClient } from "@v1/supabase/server";
import { embed, embedMany } from "ai";

const embeddingModel = openai.embedding("text-embedding-ada-002");

const generateChunks = (input: string): string[] => {
  return input
    .trim()
    .split(".")
    .filter((i) => i !== "");
};

export const generateEmbeddings = async (
  value: string,
): Promise<Array<{ embedding: number[]; content: string }>> => {
  const chunks = generateChunks(value);
  const { embeddings } = await embedMany({
    model: embeddingModel,
    values: chunks,
  });
  return embeddings.map((e, i) => ({ content: chunks[i] ?? "", embedding: e }));
};

export const generateEmbedding = async (value: string): Promise<number[]> => {
  const input = value.replaceAll("\\n", " ");
  const { embedding } = await embed({
    model: embeddingModel,
    value: input,
  });
  return embedding;
};

export const findRelavantKnowledge = async (userQuery: string) => {
  const supabase = createClient();
  const userQueryEmbedded = await generateEmbedding(userQuery);

  const { data, error } = await supabase.rpc(
    "cosine_similarity_search_knowledge",
    {
      query_embedding: userQueryEmbedded as unknown as string, //TODO: this is very bad practice,
    },
  );
  if (error) {
    console.error("Error finding similar knowledge:", error);
    return [];
  }
  return data;
};
