import { z } from "zod";

export const createMemberSchema = z.object({
  fullName: z.string().min(1),
  roleDescription: z.string().min(1),
  phoneNumber: z.string().min(1),
});
