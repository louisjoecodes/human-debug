"use server";

import { createCaseSchema } from "@/actions/cases/schema";
import { authActionClient } from "@/actions/safe-action";
import { createCase } from "@v1/supabase/mutations";
import { revalidateTag } from "next/cache";

export const createCaseAction = authActionClient
  .schema(createCaseSchema)
  .metadata({
    name: "create-case",
  })
  .action(
    async ({
      parsedInput: { first_name, last_name, date_of_birth },
      ctx: { user },
    }) => {
      const { data } = await createCase({
        first_name,
        last_name,
        date_of_birth,
        user_id: user.id,
      });
      revalidateTag("cases");
      return data;
    }
  );
