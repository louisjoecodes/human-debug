"use server";

import { createKnowledgeSchema } from "@/actions/knowledge/schema";
import { authActionClient } from "@/actions/safe-action";
import { createKnowledge } from "@v1/supabase/mutations";
import { revalidateTag } from "next/cache";

export const createKnowledgeAction = authActionClient
  .schema(createKnowledgeSchema)
  .metadata({
    name: "create-knowledge",
  })
  .action(async ({ parsedInput: { content }, ctx: { user } }) => {
    const { data } = await createKnowledge({
      content,
      user_id: user.id,
    });
    revalidateTag("knowledge");
    return data;
  });
