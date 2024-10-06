"use client";

import { useState, useEffect } from "react";
import { createClient } from "@v1/supabase/client";
import { AnalysisChat } from "./analysis-chat";
import {
    Card,
    CardContent,
    CardFooter,
    CardHeader,
    CardTitle,
} from "@v1/ui/card";
import { Heatmap } from "@/components/heatmap";
import { env } from "@/env.mjs";

export function AnalysisViewer({ caseId }: { caseId: string }) {
    const [patientData, setPatientData] = useState<any>(null);
    const [phenotypes, setPhenotypes] = useState<any[]>([]);
    const [variants, setVariants] = useState<any[]>([]);

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

    return (
        <div>
            <h1>Analysis Viewer</h1>
            <div className="grid flex-1 items-start gap-4 p-4 sm:px-6 sm:py-0 md:gap-8 lg:grid-cols-3 xl:grid-cols-3">
                <Heatmap phenotypes={phenotypes} variants={variants} />
                <div className="grid auto-rows-max items-start gap-4 md:gap-8 lg:col-span-2">
                    <div>
                        {patientData && <div>{patientData.name}</div>}
                        {phenotypes && <div>{JSON.stringify(phenotypes)}</div>}
                        {variants && <div>{JSON.stringify(variants)}</div>}
                    </div>
                </div>
                <div>
                    <Card className="overflow-hidden">
                        <CardHeader className="flex flex-row items-start bg-muted/50">
                            <div className="grid gap-0.5">
                                <CardTitle className="group flex items-center gap-2 text-lg">
                                    Analysis Chat
                                </CardTitle>
                            </div>
                        </CardHeader>
                        <CardContent className="p-6">
                            <AnalysisChat patientData={patientData} phenotypes={phenotypes} />
                        </CardContent>
                        <CardFooter className="flex flex-row items-center border-t bg-muted/50 px-6 py-3">
                            <div className="text-xs text-muted-foreground">Experimental</div>
                        </CardFooter>
                    </Card>
                </div>
            </div>
        </div>
    );
}
