import { z } from "zod";

export const shareLinkSchema = z.object({
  knowledgeId: z.string(),
  baseUrl: z.string(),
});

export const createCaseSchema = z.object({
  first_name: z.string().min(1),
  last_name: z.string().min(1),
  date_of_birth: z.string().min(1),
});
export const createReportSchema = z.object({
  caseId: z.string(),
  report: z.string(),
});

export const deleteCaseSchema = z.object({
  id: z.string(),
});
