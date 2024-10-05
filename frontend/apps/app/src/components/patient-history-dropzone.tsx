"use client";

import React, { useCallback } from 'react';
import { useDropzone } from 'react-dropzone';
import { FiUpload } from 'react-icons/fi';


export const PatientHistoryDropzone: React.FC = () => {
    const onDrop = useCallback((acceptedFiles: File[]) => {
        console.log('accepted files', acceptedFiles);
    }, []);

    const { getRootProps, getInputProps, isDragActive } = useDropzone({
        onDrop,
        accept: {
            'image/*': ['.jpeg', '.jpg', '.png'],
            'application/pdf': ['.pdf'],
        },
        multiple: true,
    });

    return (
        <div
            {...getRootProps()}
            className={`border-2 border-dashed rounded-lg p-8 text-center cursor-pointer transition-colors ${isDragActive ? 'border-blue-500 bg-blue-50' : 'border-gray-300 hover:border-gray-400'
                }`}
        >
            <input {...getInputProps()} />
            <FiUpload className="mx-auto text-4xl text-gray-400 mb-4" />
            {isDragActive ? (
                <p className="text-blue-500">Drop the files here...</p>
            ) : (
                <div>
                    <p className="text-gray-600">Drag and drop files here, or click to select files</p>
                    <p className="text-sm text-gray-400 mt-2">Supported formats: JPEG, PNG, PDF</p>
                </div>
            )}
        </div>
    );
};
