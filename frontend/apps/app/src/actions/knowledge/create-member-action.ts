"use server";

import { createMemberSchema } from "@/actions/member/schema";
import { authActionClient } from "@/actions/safe-action";
import { createMember } from "@v1/supabase/mutations";
import { revalidateTag } from "next/cache";

export const createKnowledgeAction = authActionClient
  .schema(createMemberSchema)
  .metadata({
    name: "create-member",
  })
  .action(
    async ({
      parsedInput: { fullName, roleDescription, phoneNumber },
      ctx: { user },
    }) => {
      const { data } = await createMember({
        full_name: fullName,
        role_description: roleDescription,
        phone_number: phoneNumber,
        user_id: user.id,
      });
      revalidateTag("members");
      return data;
    },
  );
