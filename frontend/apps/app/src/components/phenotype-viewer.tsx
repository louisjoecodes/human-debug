"use client"

import { useEffect, useState } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@v1/ui/card"
import { Avatar, AvatarFallback, AvatarImage } from "@v1/ui/avatar"
import { Badge } from "@v1/ui/badge"
import { ScrollArea } from "@v1/ui/scroll-area"
import { createClient } from "@v1/supabase/client"

const COLORS = [
    'hsl(var(--chart-1))',
    'hsl(var(--chart-2))',
    'hsl(var(--chart-3))',
    'hsl(var(--chart-4))',
    'hsl(var(--chart-5))',
    'hsl(var(--chart-6))',
    'hsl(var(--chart-7))',
    'hsl(var(--chart-8))',
    'hsl(var(--chart-9))',
    'hsl(var(--chart-10))',
]

export function PhenotypeViewer({ caseId }: { caseId: string }) {
    const [patientData, setPatientData] = useState<any>(null)
    const [phenotypes, setPhenotypes] = useState<any[]>([])


    useEffect(() => {
        async function fetchPatientData() {
            const supabase = createClient()
            const { data, error } = await supabase
                .from('reports')
                .select('content_json, created_at')
                .eq('case_id', caseId)
                .order('created_at', { ascending: true })

            if (error) {
                console.error('Error fetching patient data:', error)
            } else if (data && data.length > 0) {
                setPatientData(data[0].content_json)

                const uniquePhenotypes = new Map()
                data.forEach(report => {
                    report.content_json.phenotype_classes.forEach(phenotype => {
                        uniquePhenotypes.set(phenotype.id, {
                            ...phenotype,
                            created_at: report.created_at
                        })
                    })
                })

                setPhenotypes(Array.from(uniquePhenotypes.values()))
            }
        }

        fetchPatientData()
    }, [caseId])

    if (!patientData) {
        return <div>Loading...</div>
    }

    return (
        <div className="p-4">
            <Card className="mb-6">
                <CardHeader>
                    <CardTitle>Phenotype Information</CardTitle>
                </CardHeader>
                <CardContent>
                    <ScrollArea className="h-[700px] w-full">
                        {phenotypes.map((phenotype, index) => (
                            <div key={phenotype.id} className="mb-4 p-4 border rounded-lg">
                                <h3 className="text-lg font-semibold flex items-center">
                                    <span className="w-4 h-4 rounded-full mr-2" style={{ backgroundColor: COLORS[index % COLORS.length] }}></span>
                                    {phenotype.name}
                                </h3>
                                <p className="text-sm text-muted-foreground mt-1">{phenotype.definition}</p>
                                <p className="text-xs text-muted-foreground mt-1">ID: {phenotype.id}</p>
                                <p className="text-xs text-muted-foreground mt-1">Created: {new Date(phenotype.created_at).toLocaleString()}</p>
                            </div>
                        ))}
                    </ScrollArea>
                </CardContent>
            </Card>
            <Card>
                <CardHeader>
                    <CardTitle>Patient Profile</CardTitle>
                </CardHeader>
                <CardContent className="flex items-center space-x-4">
                    <Avatar className="h-24 w-24">
                        <AvatarImage src="/placeholder.svg?height=96&width=96" alt={`${patientData.patient.first_name} ${patientData.patient.last_name}`} />
                        <AvatarFallback>{patientData.patient.first_name[0]}{patientData.patient.last_name[0]}</AvatarFallback>
                    </Avatar>
                    <div>
                        <h2 className="text-2xl font-semibold">{patientData.patient.first_name} {patientData.patient.last_name}</h2>
                        <p className="text-muted-foreground">Age: {patientData.patient.age}</p>
                        <p className="text-muted-foreground">Gender: {patientData.patient.gender}</p>
                        <p className="text-muted-foreground">DOB: {patientData.patient.date_of_birth}</p>
                        <Badge variant="destructive" className="mt-2">{patientData.patient.disease}</Badge>
                    </div>
                </CardContent>
            </Card>



        </div>
    )
}