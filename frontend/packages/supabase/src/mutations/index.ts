import { logger } from "@v1/logger";
import { createClient } from "@v1/supabase/server";
import { generateEmbeddings } from "../lib/ai/embedding";
import type { Database, Tables, TablesInsert, TablesUpdate } from "../types";

export async function updateUser(userId: string, data: TablesUpdate<"users">) {
  const supabase = createClient();

  try {
    const result = await supabase.from("users").update(data).eq("id", userId);

    return result;
  } catch (error) {
    logger.error(error);

    throw error;
  }
}

export async function createCase(caseInsert: TablesInsert<"cases">) {
  const supabase = createClient();

  console.log("CRETING CASE");
  try {
    // Insert Case
    const { data: createdCase, error: insertError } = await supabase
      .from("cases")
      .insert(caseInsert)
      .select()
      .single();

    if (insertError) throw insertError;
    console.log(createdCase);
    console.log(insertError);

    // // Generate embeddings
    // const embeddings = await generateEmbeddings(insertedKnowledge.content);

    // // Insert embeddings
    // const { error: embeddingError } = await supabase
    //   .from("knowledge_content_embeddings")
    //   .insert(
    //     embeddings.map((embedding) => ({
    //       knowledge_id: insertedKnowledge.id,
    //       content: embedding.content,
    //       embedding: JSON.stringify(embedding.embedding),
    //     }))
    //   );
    // console.log(embeddingError);

    // if (embeddingError) throw embeddingError;

    return {
      message: "Case successfully created.",
      data: createdCase,
    };
  } catch (error) {
    logger.error(error);
    return { message: "Error creating knowledge and embeddings", error };
  }
}

export const uploadReport = async (file: File) => {
  const fileName = `${Date.now()}_${file.name}`;
  const { data, error } = await supabase.storage
    .from("reports")
    .upload(fileName, file);

  return { path: data?.path, error };
};

export const createReport = async ({
  case_id,
  file_path,
  created_by,
}: {
  case_id: string;
  file_path: string;
  created_by: string;
}) => {
  return await supabase
    .from("reports")
    .insert({
      case_id,
      file_path,
      created_by,
    })
    .select()
    .single();
};
