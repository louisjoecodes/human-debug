"use server";

import { createMemberSchema } from "@/actions/member/schema";
import { authActionClient } from "@/actions/safe-action";
import { createMember } from "@v1/supabase/mutations";
import { revalidateTag } from "next/cache";

export const createMemberAction = authActionClient
  .schema(createMemberSchema)
  .metadata({
    name: "create-user",
  })
  .action(async ({ parsedInput }) => {
    const { data } = await createMember(parsedInput);
    revalidateTag("users");
    return data;
  });
