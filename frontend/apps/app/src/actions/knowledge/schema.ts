import { z } from "zod";

export const shareLinkSchema = z.object({
  knowledgeId: z.string(),
  baseUrl: z.string(),
});

export const createKnowledgeSchema = z.object({
  content: z.string().min(1),
});

export const deleteKnowledgeSchema = z.object({
  id: z.string(),
});
