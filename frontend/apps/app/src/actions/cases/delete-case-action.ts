"use server";

import { deleteCaseSchema } from "@/actions/cases/schema";
import { authActionClient } from "@/actions/safe-action";
import { createCase } from "@v1/supabase/mutations";
import { createClient } from "@v1/supabase/server";
import { revalidateTag } from "next/cache";

export const deleteCaseAction = authActionClient
  .schema(deleteCaseSchema)
  .metadata({
    name: "delete-case",
  })
  .action(async ({ parsedInput: { id }, ctx: { user } }) => {
    const supabase = await createClient();
    await supabase.from("cases").delete().eq("id", id);
    revalidateTag("cases");
    return id;
  });
