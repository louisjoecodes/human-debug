"use client";

import { useState, useEffect } from "react";
import { createClient } from "@v1/supabase/client";
import { AnalysisChat } from "./analysis-chat";
import {
    Card,
    CardContent,
    CardHeader,
    CardTitle,
} from "@v1/ui/card";
import { env } from "@/env.mjs";
import ReactMarkdown from 'react-markdown';
import { Tabs, TabsList, TabsTrigger, TabsContent } from "@v1/ui/tabs";
import { Heatmap } from "./heatmap";
import { Button } from "@v1/ui/button";
import { Loader2 } from "lucide-react";

export function AnalysisViewer({ caseId }: { caseId: string }) {
    const [patientData, setPatientData] = useState<any>(null);
    const [phenotypes, setPhenotypes] = useState<any[]>([]);
    const [variants, setVariants] = useState<any[]>([]);
    const [activeTab, setActiveTab] = useState<'search' | 'analysis'>('search');
    const [isGeneratingReport, setIsGeneratingReport] = useState(false);
    const [reportContent, setReportContent] = useState<string>('');

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

    const generateReport = async () => {
        setIsGeneratingReport(true);

        // Simulate API call delay
        await new Promise(resolve => setTimeout(resolve, 5000));

        const newReportContent = `
# Patient Report: Maria Isabella Rodriguez

Patient ID: MIR2024001

## Clinical Presentation

Maria Isabella Rodriguez, a 32-year-old female, presents with a complex array of symptoms suggesting a potential genetic disorder. The patient reports a history of progressive muscle weakness, recurrent respiratory infections, and delayed motor development since childhood.

## Phenotypic Analysis

- Muscular hypotonia
- Recurrent lower respiratory tract infections
- Delayed gross motor development
- Mild intellectual disability
- Distinctive facial features: high-arched palate, long face

## Genetic Analysis

Whole Exome Sequencing (WES) was performed, revealing the following significant variant:

- Gene: DMD (Dystrophin)
- Variant: c.6638G>T (p.Gly2213Val)
- Zygosity: Heterozygous
- Classification: Likely Pathogenic

## Interpretation

The identified variant in the DMD gene is consistent with a diagnosis of Duchenne Muscular Dystrophy (DMD) or Becker Muscular Dystrophy (BMD). The heterozygous state in this female patient suggests she may be a carrier, but the presence of symptoms indicates a potential case of manifesting carrier syndrome.

## Recommendations

- Referral to a neuromuscular specialist for comprehensive evaluation
- Cardiac assessment to monitor for potential cardiomyopathy
- Genetic counseling for the patient and family members
- Consider muscle biopsy for dystrophin protein analysis
- Develop a multidisciplinary care plan including physical therapy, respiratory support, and regular monitoring

This report was generated based on the current available information and should be interpreted in conjunction with clinical findings and additional diagnostic tests.
        `;

        setReportContent(newReportContent);
        setIsGeneratingReport(false);
    };

    const MarkdownComponents = {
        h1: ({ node, ...props }) => <h1 className="text-3xl font-bold mt-6 mb-4" {...props} />,
        h2: ({ node, ...props }) => <h2 className="text-2xl font-semibold mt-5 mb-3" {...props} />,
        h3: ({ node, ...props }) => <h3 className="text-xl font-medium mt-4 mb-2" {...props} />,
        p: ({ node, ...props }) => <p className="mb-4" {...props} />,
        ul: ({ node, ...props }) => <ul className="list-disc pl-6 mb-4" {...props} />,
        ol: ({ node, ...props }) => <ol className="list-decimal pl-6 mb-4" {...props} />,
        li: ({ node, ...props }) => <li className="mb-1" {...props} />,
        strong: ({ node, ...props }) => <strong className="font-semibold" {...props} />,
        blockquote: ({ node, ...props }) => (
            <blockquote className="border-l-4 border-gray-300 pl-4 italic my-4" {...props} />
        ),
    };

    return (
        <div className="flex flex-col h-[calc(100vh-24px)]"> {/* Adjust 64px if your header height is different */}
            <div className="grid flex-1 items-start gap-4 p-4 sm:px-6 sm:py-4 md:gap-8 lg:grid-cols-3 xl:grid-cols-3 h-full overflow-hidden">
                <Card className="lg:col-span-2 xl:col-span-2 flex flex-col h-full max-h-[calc(100vh-100px)]">
                    <CardHeader className="bg-muted/50">
                        <CardTitle className="text-lg">
                            <div className="flex justify-between items-center">
                                Patient Case Report
                                <Button
                                    variant="primary"
                                    size="sm"
                                    onClick={generateReport}
                                    disabled={isGeneratingReport}
                                >
                                    {isGeneratingReport ? (
                                        <>
                                            <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                                            Generating...
                                        </>
                                    ) : (
                                        'Generate Report 🧬'
                                    )}
                                </Button>
                            </div>
                        </CardTitle>
                    </CardHeader>
                    <CardContent className="p-4 flex-1 overflow-y-auto">
                        <ReactMarkdown
                            components={MarkdownComponents}
                            className="prose dark:prose-invert max-w-none"
                        >
                            {reportContent}
                        </ReactMarkdown>
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
                <Card className="flex flex-col col-span-3 h-full">
                    <CardHeader>Phenotype vs. Genotype Correlation</CardHeader>
                    <CardContent>
                        <Heatmap />
                    </CardContent>
                </Card>

            </div>


        </div>
    );
}
