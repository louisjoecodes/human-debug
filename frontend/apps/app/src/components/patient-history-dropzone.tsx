"use client";

import React, { useCallback, useState } from 'react';
import { useDropzone } from 'react-dropzone';
import { FiUpload } from 'react-icons/fi';
import axios from 'axios';
import { createReportAction } from '@/actions/cases/reports/create-report-action';
import { createClient } from "@v1/supabase/client";
import { nanoid } from 'nanoid';

export const PatientHistoryDropzone = ({ caseId }: { caseId: string }) => {
    const [patientInfo, setPatientInfo] = useState(null);
    const [isLoading, setIsLoading] = useState(false);

    const onDrop = useCallback(async (acceptedFiles: File[]) => {
        setIsLoading(true);
        const file = acceptedFiles[0];

        try {
            const supabase = createClient();
            const id = nanoid()
            if (!file) {
                throw new Error('No file selected');
            }

            const { data, error } = await supabase.storage
                .from('reports')
                .upload(`${caseId}/${id}`, file);

            const formData = new FormData();
            formData.append('file', file);

            const response = await axios.post('http://localhost:8000/process_letter', formData, {
                headers: {
                    'Content-Type': 'multipart/form-data',
                },
            });
            setPatientInfo(response.data);

            // Call createReportAction to add the file to the database
            const result = await createReportAction({
                caseId,
                report: file,
            });
            if (!result) {
                throw new Error('Failed to create report');
            }

            if (result.data) {
                console.log('Report created successfully:', result.data);
            }
        } catch (error) {
            console.error('Error processing file:', error);
        } finally {
            setIsLoading(false);
        }
    }, [caseId]);

    const { getRootProps, getInputProps, isDragActive } = useDropzone({
        onDrop,
        accept: {
            'image/*': ['.jpeg', '.jpg', '.png'],
            'application/pdf': ['.pdf'],
        },
        multiple: false,
    });

    return (
        <div>
            <div
                {...getRootProps()}
                className={`border-2 border-dashed rounded-lg p-8 text-center cursor-pointer transition-colors ${isDragActive ? 'border-blue-500 bg-blue-50' : 'border-gray-300 hover:border-gray-400'
                    }`}
            >
                <input {...getInputProps()} />
                <FiUpload className="mx-auto text-4xl text-gray-400 mb-4" />
                {isDragActive ? (
                    <p className="text-blue-500">Drop the file here...</p>
                ) : (
                    <div>
                        <p className="text-gray-600">Drag and drop a file here, or click to select a file</p>
                        <p className="text-sm text-gray-400 mt-2">Supported formats: JPEG, PNG, PDF</p>
                    </div>
                )}
            </div>
            {isLoading && <p className="mt-4 text-center">Processing file...</p>}
            {patientInfo && (
                <div className="mt-6 p-4 border rounded-lg">
                    <h3 className="text-lg font-semibold mb-2">Extracted Patient Information:</h3>
                    <pre className="whitespace-pre-wrap">{JSON.stringify(patientInfo, null, 2)}</pre>
                </div>
            )}
        </div>
    );
};
