"use server";

import { deleteKnowledgeSchema } from "@/actions/knowledge/schema";
import { authActionClient } from "@/actions/safe-action";
import { createKnowledge } from "@v1/supabase/mutations";
import { createClient } from "@v1/supabase/server";
import { revalidateTag } from "next/cache";

export const deleteKnowledgeAction = authActionClient
  .schema(deleteKnowledgeSchema)
  .metadata({
    name: "delete-knowledge",
  })
  .action(async ({ parsedInput: { id }, ctx: { user } }) => {
    const supabase = await createClient();
    await supabase.from("knowledge").delete().eq("id", id);
    revalidateTag("knowledge");
    return id;
  });
