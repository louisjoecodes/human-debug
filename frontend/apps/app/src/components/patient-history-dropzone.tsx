"use client";

import React, { useCallback, useState } from "react";
import { useDropzone } from "react-dropzone";
import { FiUpload } from "react-icons/fi";
import axios from "axios";
import { createReportAction } from "@/actions/cases/reports/create-report-action";
import { createClient } from "@v1/supabase/client";
import { nanoid } from "nanoid";
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogFooter,
  DialogClose,
} from "@v1/ui/dialog";
import { Button } from "@v1/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@v1/ui/card";
import { ScrollArea } from "@v1/ui/scroll-area";
import { env } from "@/env.mjs";

const COLORS = [
  "hsl(var(--chart-1))",
  "hsl(var(--chart-2))",
  "hsl(var(--chart-3))",
  "hsl(var(--chart-4))",
  "hsl(var(--chart-5))",
  "hsl(var(--chart-6))",
  "hsl(var(--chart-7))",
  "hsl(var(--chart-8))",
  "hsl(var(--chart-9))",
  "hsl(var(--chart-10))",
];

export const PatientHistoryDropzone = ({ caseId }: { caseId: string }) => {
  const [patientInfo, setPatientInfo] = useState(null);
  const [isLoading, setIsLoading] = useState(false);
  const [isDialogOpen, setIsDialogOpen] = useState(false);
  const [ocrContent, setOcrContent] = useState("");
  const [phenotypes, setPhenotypes] = useState<any[]>([]);

  const onDrop = useCallback(
    async (acceptedFiles: File[]) => {
      setIsLoading(true);
      const file = acceptedFiles[0];

      try {
        const supabase = createClient();
        const id = nanoid();
        if (!file) {
          throw new Error("No file selected");
        }

        const { data, error } = await supabase.storage
          .from("reports")
          .upload(`${caseId}/${id}`, file);

        const formData = new FormData();
        formData.append("file", file);

        const response = await axios.post(
          `${env.NEXT_PUBLIC_HUMAN_DEBUG_BACKEND_API_URL}/process_letter`,
          formData,
          {
            headers: {
              "Content-Type": "multipart/form-data",
            },
          },
        );

        setOcrContent(response.data.content);
        setIsDialogOpen(true);
        setPatientInfo(response.data);
        setPhenotypes(response.data.phenotype_classes || []);

        // Call createReportAction to add the file to the database
        const result = await createReportAction({
          caseId,
          report: file,
        });
        if (!result) {
          throw new Error("Failed to create report");
        }

        if (result.data) {
          console.log("Report created successfully:", result.data);
        }
      } catch (error) {
        console.error("Error processing file:", error);
      } finally {
        setIsLoading(false);
      }
    },
    [caseId],
  );

  const handleConfirmOCR = async () => {
    try {
      setIsLoading(true);
      const supabase = createClient();
      const { data: insertedData, error } = await supabase
        .from("reports")
        .insert({
          case_id: caseId,
          content_json: patientInfo,
        });

      if (error) throw error;

      console.log("Report saved successfully:", insertedData);
      setIsDialogOpen(false);
      Router.refresh();
    } catch (error) {
      console.error("Error saving report:", error);
    } finally {
      setIsLoading(false);
    }
  };

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    accept: {
      "image/*": [".jpeg", ".jpg", ".png"],
      "application/pdf": [".pdf"],
    },
    multiple: false,
  });

  return (
    <div>
      <div
        {...getRootProps()}
        className={`border-2 border-dashed rounded-lg p-8 text-center cursor-pointer transition-colors ${isDragActive
          ? "border-blue-500 bg-blue-50"
          : "border-gray-300 hover:border-gray-400"
          }`}
      >
        <input {...getInputProps()} />
        <FiUpload className="mx-auto text-4xl text-gray-400 mb-4" />
        {isDragActive ? (
          <p className="text-blue-500">Drop the file here...</p>
        ) : (
          <div>
            <p className="text-gray-600">
              Drag and drop a file here, or click to select a file
            </p>
            <p className="text-sm text-gray-400 mt-2">
              Supported formats: JPEG, PNG, PDF
            </p>
          </div>
        )}
      </div>
      {isLoading && <p className="mt-4 text-center">Processing file...</p>}

      <Dialog open={isDialogOpen} onOpenChange={setIsDialogOpen}>
        <DialogContent className="max-w-4xl">
          <DialogHeader>
            <DialogTitle>Confirm patient Phenotype changes</DialogTitle>
          </DialogHeader>
          <div className="grid grid-cols-1 gap-4">
            <div>
              <Card>
                <CardHeader>
                  <CardTitle>Detected Phenotype Classifications 🧬</CardTitle>
                </CardHeader>
                <CardContent>
                  <ScrollArea className="h-[60vh]">
                    {phenotypes.map((phenotype, index) => (
                      <div
                        key={phenotype.id}
                        className="mb-4 p-4 border rounded-lg"
                      >
                        <h3 className="text-lg font-semibold flex items-center">
                          <span
                            className="w-4 h-4 rounded-full mr-2"
                            style={{
                              backgroundColor: COLORS[index % COLORS.length],
                            }}
                          ></span>
                          {phenotype.name}
                        </h3>
                        <p className="text-sm text-muted-foreground mt-1">
                          {phenotype.definition}
                        </p>
                        <p className="text-xs text-muted-foreground mt-1">
                          ID: {phenotype.id}
                        </p>
                      </div>
                    ))}
                  </ScrollArea>
                </CardContent>
              </Card>
            </div>
          </div>
          <DialogFooter>
            <Button onClick={handleConfirmOCR} disabled={isLoading}>
              {isLoading ? "Saving..." : "Confirm and Save"}
            </Button>
            <DialogClose asChild>
              <Button variant="outline" disabled={isLoading}>
                Cancel
              </Button>
            </DialogClose>
          </DialogFooter>
        </DialogContent>
      </Dialog>
    </div>
  );
};
