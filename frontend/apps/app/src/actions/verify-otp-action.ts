"use server";

import { actionClient } from "@/actions/safe-action";
import { verifyOtpSchema } from "@/actions/schema";
import { createClient } from "@v1/supabase/server";
import { redirect } from "next/navigation";

export const verifyOtpAction = actionClient
  .schema(verifyOtpSchema)
  .action(async ({ parsedInput: { email, token } }) => {
    const supabase = createClient();
    try {
      const res = await supabase.auth.verifyOtp({
        email,
        token,
        type: "email",
      });
      if (res.error) {
        throw new Error(res.error.message);
      }
    } catch (error) {
      throw new Error(error instanceof Error ? error.message : "Unknown error");
    }
    redirect("/");
  });
