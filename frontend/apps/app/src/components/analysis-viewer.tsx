"use client";

import { useState, useEffect, useRef } from "react";
import { createClient } from "@v1/supabase/client";
import { AnalysisChat } from "./analysis-chat";
import {
    Card,
    CardContent,
    CardFooter,
    CardHeader,
    CardTitle,
} from "@v1/ui/card";
import { env } from "@/env.mjs";
import dynamic from 'next/dynamic';
import { Tabs, TabsList, TabsTrigger, TabsContent } from "@v1/ui/tabs";

const EditorJS = dynamic(() => import('./EditorJS'), { ssr: false });

export function AnalysisViewer({ caseId }: { caseId: string }) {
    const [patientData, setPatientData] = useState<any>(null);
    const [phenotypes, setPhenotypes] = useState<any[]>([]);
    const [variants, setVariants] = useState<any[]>([]);
    const [activeTab, setActiveTab] = useState<'search' | 'analysis'>('search');

    useEffect(() => {
        const fetchVariants = async () => {
            try {
                const response = await fetch(`${env.NEXT_PUBLIC_HUMAN_DEBUG_BACKEND_API_URL}/variants`);
                const data = await response.json();
                setVariants(data);
            } catch (error) {
                console.error("Error fetching variants:", error);
            }
        };

        fetchVariants();
    }, []);

    useEffect(() => {
        async function fetchPatientData() {
            const supabase = createClient();
            const { data, error } = await supabase
                .from("reports")
                .select("content_json, created_at")
                .eq("case_id", caseId)
                .order("created_at", { ascending: true });

            if (error) {
                console.error("Error fetching patient data:", error);
            } else if (data && data.length > 0) {
                setPatientData(data[0].content_json);

                const uniquePhenotypes = new Map();
                data.forEach((report) => {
                    report.content_json.phenotype_classes.forEach((phenotype) => {
                        uniquePhenotypes.set(phenotype.id, {
                            ...phenotype,
                            created_at: report.created_at,
                        });
                    });
                });

                setPhenotypes(Array.from(uniquePhenotypes.values()));
            }
        }

        fetchPatientData();
    }, [caseId]);

    const defaultEditorContent = {
        time: new Date().getTime(),
        blocks: [
            {
                type: "header",
                data: {
                    text: "Patient Case Report",
                    level: 1
                }
            },
            {
                type: "paragraph",
                data: {
                    text: "Patient presents with the following phenotypes:"
                }
            },
            {
                type: "list",
                data: {
                    style: "unordered",
                    items: [
                        "Phenotype 1",
                        "Phenotype 2",
                        "Phenotype 3"
                    ]
                }
            },
            {
                type: "header",
                data: {
                    text: "Genome Sequence Analysis",
                    level: 2
                }
            },
            {
                type: "paragraph",
                data: {
                    text: "The genome sequence analysis revealed the following variants:"
                }
            },
            {
                type: "list",
                data: {
                    style: "ordered",
                    items: [
                        "Variant 1: Description",
                        "Variant 2: Description",
                        "Variant 3: Description"
                    ]
                }
            }
        ]
    };

    return (
        <div className="flex flex-col h-[calc(100vh-204px)]"> {/* Adjust 64px if your header height is different */}
            <div className="grid flex-1 items-start gap-4 p-4 sm:px-6 sm:py-4 md:gap-8 lg:grid-cols-3 xl:grid-cols-3 h-full overflow-hidden">
                <Card className="lg:col-span-2 xl:col-span-2 flex flex-col h-full">
                    <CardHeader className="bg-muted/50">
                        <CardTitle className="text-lg">Patient Case Report</CardTitle>
                    </CardHeader>
                    <CardContent className="p-0 flex-1 overflow-y-auto">
                        <EditorJS defaultValue={defaultEditorContent} />
                    </CardContent>
                </Card>

                <Card className="flex flex-col h-full">
                    <CardHeader className="bg-muted/50">
                        <CardTitle className="text-lg">Analysis Helpers</CardTitle>
                    </CardHeader>
                    <CardContent className="p-4 flex-1 overflow-hidden">
                        <Tabs value={activeTab} onValueChange={(value) => setActiveTab(value as 'search' | 'analysis')} className="flex flex-col h-full">
                            <TabsList className="grid w-full grid-cols-2">
                                <TabsTrigger value="search">Human Debug Search</TabsTrigger>
                                <TabsTrigger value="analysis">Medical Information Assistant</TabsTrigger>
                            </TabsList>
                            <TabsContent value="search" className="flex-1 overflow-hidden">
                                <div className="h-full">
                                    <iframe
                                        src="https://human-debug-search.vercel.app"
                                        title="Research Assistant"
                                        className="w-full h-full border-0"
                                    />
                                </div>
                            </TabsContent>
                            <TabsContent value="analysis" className="flex-1 overflow-hidden">
                                <AnalysisChat patientData={patientData} phenotypes={phenotypes} />
                            </TabsContent>
                        </Tabs>
                    </CardContent>
                </Card>
            </div>
        </div>
    );
}
