"use server";

import { createReportSchema } from "@/actions/cases/reports/schema";
import { authActionClient } from "@/actions/safe-action";
import { uploadReport, createReport } from "@v1/supabase/mutations";
import { revalidateTag } from "next/cache";

export const createReportAction = authActionClient
  .schema(createReportSchema)
  .metadata({
    name: "create-report",
  })
  .action(async ({ parsedInput: { caseId, report }, ctx: { user } }) => {
    try {
      // Upload the report file to the 'reports' bucket
      const { path, error: uploadError } = await uploadReport(report);

      if (uploadError) {
        throw new Error(`Failed to upload report: ${uploadError.message}`);
      }

      // Create a new report entry in the database
      const { data: reportData, error: createError } = await createReport({
        case_id: caseId,
        file_path: path,
        created_by: user.id,
      });

      if (createError) {
        throw new Error(
          `Failed to create report entry: ${createError.message}`,
        );
      }

      revalidateTag("reports");
      return reportData;
    } catch (error) {
      console.error("Error in createReportAction:", error);
      throw error;
    }
  });
