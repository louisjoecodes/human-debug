"use client"

import { useState, useEffect } from "react"
import { Button } from "@v1/ui/button"
import { ChevronRight } from "lucide-react"
import dynamic from 'next/dynamic'
import { createClient } from "@v1/supabase/client"

const Dialog = dynamic(() => import("@v1/ui/dialog").then(mod => mod.Dialog), { ssr: false })
const DialogTrigger = dynamic(() => import("@v1/ui/dialog").then(mod => mod.DialogTrigger), { ssr: false })
const DialogContent = dynamic(() => import("@v1/ui/dialog").then(mod => mod.DialogContent), { ssr: false })
const DialogHeader = dynamic(() => import("@v1/ui/dialog").then(mod => mod.DialogHeader), { ssr: false })
const DialogTitle = dynamic(() => import("@v1/ui/dialog").then(mod => mod.DialogTitle), { ssr: false })

export const ReportViewer = ({ caseId, data }: { caseId: string, data: any }) => {
    const [pdfUrl, setPdfUrl] = useState<string | null>(null)

    useEffect(() => {
        const fetchSignedUrl = async () => {
            if (caseId && data && data.name) {
                const supabase = createClient()
                const { data: signedUrl, error } = await supabase.storage
                    .from('reports')
                    .createSignedUrl(`${caseId}/${data.name}`, 60 * 60) // 1 hour expiration

                if (error) {
                    console.error('Error creating signed URL:', error)
                    return
                }

                setPdfUrl(signedUrl.signedUrl)
            }
        }

        fetchSignedUrl()
    }, [caseId, data])

    return (
        <Dialog>
            <DialogTrigger asChild>
                <Button variant="ghost" size="sm">
                    View Report
                    <ChevronRight className="ml-2 h-4 w-4" />
                </Button>
            </DialogTrigger>
            <DialogContent className="max-w-4xl h-[90vh] p-0">
                <DialogHeader className="p-4">
                    <DialogTitle>Report Viewer</DialogTitle>
                </DialogHeader>
                <div className="h-[calc(90vh-4rem)] overflow-hidden">
                    {pdfUrl ? (
                        <iframe
                            src={`${pdfUrl}#toolbar=0`}
                            className="w-full h-full"
                            title="Report PDF"
                        />
                    ) : (
                        <div className="flex items-center justify-center h-full">Loading report...</div>
                    )}
                </div>
            </DialogContent>
        </Dialog>
    )
}