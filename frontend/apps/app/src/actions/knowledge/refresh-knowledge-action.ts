"use server";

import { authActionClient } from "@/actions/safe-action";
import { revalidateTag } from "next/cache";
import { z } from "zod";

export const refreshKnowledgeAction = authActionClient
  .schema(z.object({}))
  .metadata({
    name: "refresh-knowledge",
  })
  .action(async () => {
    revalidateTag("knowledge");
  });
