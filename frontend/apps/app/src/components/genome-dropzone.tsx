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
import { Progress } from "@v1/ui/progress";

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

interface GenomeDropzoneProps {
    setProcessingComplete: (value: boolean) => void;
}

export const GenomeDropzone: React.FC<GenomeDropzoneProps> = ({ setProcessingComplete }) => {
    const [isLoading, setIsLoading] = useState(false);
    const [progress, setProgress] = useState(0);
    const [currentStep, setCurrentStep] = useState("");

    const simulateProcessing = useCallback(async () => {
        const steps = [
            "Reading FASTQ file",
            "Aligning sequences",
            "Identifying variants",
            "Generating report",
        ];

        for (let i = 0; i < steps.length; i++) {
            setCurrentStep(steps[i]);
            const stepProgress = Math.random() * 15 + 10; // Random progress between 10-25 per step
            for (let j = 0; j < stepProgress; j++) {
                await new Promise(resolve => setTimeout(resolve, 40));
                setProgress(prev => Math.min(prev + 1, 100));
            }
        }

        // Ensure we reach 100% and wait for a total of 4 seconds
        setProgress(100);
        await new Promise(resolve => setTimeout(resolve, 4000 - (100 * 40)));
        setProcessingComplete(true);
        setIsLoading(false);
    }, [setProcessingComplete]);

    const onDrop = useCallback(
        async (acceptedFiles: File[]) => {
            const file = acceptedFiles[0];
            if (file.name.endsWith('.fastq')) {
                setIsLoading(true);
                setProgress(0);
                simulateProcessing();
            } else {
                alert('Please upload a .fastq file');
            }
        },
        [simulateProcessing]
    );

    const { getRootProps, getInputProps, isDragActive } = useDropzone({
        onDrop,
        accept: {
            "application/octet-stream": [".fastq"],
        },
        multiple: false,
    });

    return (
        <div>
            {!isLoading ? (
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
                        <p className="text-blue-500">Drop the .fastq file here...</p>
                    ) : (
                        <div>
                            <p className="text-gray-600">
                                Drag and drop a .fastq file here, or click to select a file
                            </p>
                            <p className="text-sm text-gray-400 mt-2">
                                Supported format: .fastq
                            </p>
                        </div>
                    )}
                </div>
            ) : (
                <div className="mt-4">
                    <p className="text-center mb-2 font-semibold">{currentStep}</p>
                    <Progress value={progress} className="w-full h-2 bg-green-100">
                        <div
                            className="h-full bg-green-500 transition-all duration-300 ease-out"
                            style={{ width: `${progress}%` }}
                        />
                    </Progress>
                </div>
            )}
        </div>
    );
};
