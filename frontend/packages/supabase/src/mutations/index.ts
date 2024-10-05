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

export async function createKnowledge(knowledge: TablesInsert<"knowledge">) {
  const supabase = createClient();

  try {
    // Insert knowledge
    const { data: insertedKnowledge, error: insertError } = await supabase
      .from("knowledge")
      .insert(knowledge)
      .select()
      .single();

    if (insertError) throw insertError;

    // Generate embeddings
    const embeddings = await generateEmbeddings(insertedKnowledge.content);

    // Insert embeddings
    const { error: embeddingError } = await supabase
      .from("knowledge_content_embeddings")
      .insert(
        embeddings.map((embedding) => ({
          knowledge_id: insertedKnowledge.id,
          content: embedding.content,
          embedding: JSON.stringify(embedding.embedding),
        })),
      );
    console.log(embeddingError);

    if (embeddingError) throw embeddingError;

    return {
      message: "Knowledge successfully created and embedded.",
      data: insertedKnowledge,
    };
  } catch (error) {
    logger.error(error);
    return { message: "Error creating knowledge and embeddings", error };
  }
}

export async function createMember(member: TablesInsert<"members">) {
  const supabase = createClient();

  try {
    // Insert member
    const { data: insertedMember, error: insertError } = await supabase
      .from("members")
      .insert(member)
      .select()
      .single();

    if (insertError) throw insertError;

    // Generate embeddings
    const embeddings = await generateEmbeddings(
      insertedMember.role_description,
    );

    // Insert embeddings
    const { error: embeddingError } = await supabase
      .from("members_role_description_embeddings")
      .insert(
        embeddings.map(({ embedding, content }) => ({
          member_id: insertedMember.id,
          role_description: content,
          embedding: JSON.stringify(embedding),
        })),
      );
    console.log(embeddingError);

    if (embeddingError) throw embeddingError;

    return {
      message: "Member successfully created and embedded.",
      data: insertedMember,
    };
  } catch (error) {
    logger.error(error);
    return { message: "Error creating member and embeddings", error };
  }
}
