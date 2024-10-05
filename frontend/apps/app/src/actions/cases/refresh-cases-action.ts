"use server";

import { authActionClient } from "@/actions/safe-action";
import { revalidateTag } from "next/cache";
import { z } from "zod";

export const refreshCasesAction = authActionClient
  .schema(z.object({}))
  .metadata({
    name: "refresh-cases",
  })
  .action(async () => {
    revalidateTag("cases");
  });
