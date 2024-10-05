import { z } from "zod";

export const createReportSchema = z.object({
  caseId: z.string(),
  report: z.instanceof(File),
});
