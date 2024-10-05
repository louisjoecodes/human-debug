import { getCaseById } from "@v1/supabase/queries";
import { Avatar, AvatarFallback, AvatarImage } from "@v1/ui/avatar";
import { Button } from "@v1/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@v1/ui/card";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@v1/ui/table";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@v1/ui/tabs";
import { ChevronRight, Stethoscope } from "lucide-react";
import { PatientHistoryDropzone } from "@/components/patient-history-dropzone";

export default async function Page({ params }: { params: { id: string } }) {
    const { data } = await getCaseById(params.id);
    if (!data) {
        return <div>Case not found</div>;
    }

    const { first_name, last_name, date_of_birth, id } = data;

    return (
        <>
            <header className="mb-8">
                <h1 className="text-3xl font-bold">{first_name} {last_name}</h1>
                <p className="text-gray-500">Date of Birth: {date_of_birth}</p>
            </header>
            <Tabs defaultValue="overview">
                <TabsList className="mb-6">
                    <TabsTrigger value="overview">Overview</TabsTrigger>
                    <TabsTrigger value="appointments">Appointments</TabsTrigger>
                    <TabsTrigger value="medications">Medications</TabsTrigger>
                    <TabsTrigger value="notes">Notes</TabsTrigger>
                </TabsList>
                <TabsContent value="overview">
                    <Card className="mb-6">
                        <CardHeader>
                            <CardTitle>Upload Patient History</CardTitle>
                        </CardHeader>
                        <CardContent>
                            <PatientHistoryDropzone caseId={id} />
                        </CardContent>
                    </Card>
                    <Card className="md:col-span-2">
                        <CardHeader>
                            <CardTitle>Recent Medical History</CardTitle>
                        </CardHeader>
                        <CardContent>
                            <Table>
                                <TableHeader>
                                    <TableRow>
                                        <TableHead>Date</TableHead>
                                        <TableHead>Description</TableHead>
                                        <TableHead>Doctor</TableHead>
                                        <TableHead>Action</TableHead>
                                    </TableRow>
                                </TableHeader>
                                <TableBody>
                                    <TableRow>
                                        <TableCell>05/15/2023</TableCell>
                                        <TableCell>Annual physical examination</TableCell>
                                        <TableCell>Dr. Sarah Johnson</TableCell>
                                        <TableCell>
                                            <Button variant="ghost" size="sm">
                                                View Details
                                                <ChevronRight className="ml-2 h-4 w-4" />
                                            </Button>
                                        </TableCell>
                                    </TableRow>
                                    <TableRow>
                                        <TableCell>03/22/2023</TableCell>
                                        <TableCell>Flu vaccination</TableCell>
                                        <TableCell>Dr. Michael Lee</TableCell>
                                        <TableCell>
                                            <Button variant="ghost" size="sm">
                                                View Details
                                                <ChevronRight className="ml-2 h-4 w-4" />
                                            </Button>
                                        </TableCell>
                                    </TableRow>
                                    <TableRow>
                                        <TableCell>01/10/2023</TableCell>
                                        <TableCell>Sprained ankle treatment</TableCell>
                                        <TableCell>Dr. Emily Chen</TableCell>
                                        <TableCell>
                                            <Button variant="ghost" size="sm">
                                                View Details
                                                <ChevronRight className="ml-2 h-4 w-4" />
                                            </Button>
                                        </TableCell>
                                    </TableRow>
                                </TableBody>
                            </Table>
                        </CardContent>
                    </Card>

                </TabsContent>
                <TabsContent value="appointments">
                    <Card>
                        <CardHeader>
                            <CardTitle>Upcoming Appointments</CardTitle>
                        </CardHeader>
                        <CardContent>
                            <div className="space-y-6">
                                <div className="flex items-center space-x-4">
                                    <Avatar>
                                        <AvatarImage src="/placeholder.svg?height=40&width=40" alt="Dr. Robert Smith" />
                                        <AvatarFallback>RS</AvatarFallback>
                                    </Avatar>
                                    <div className="flex-1">
                                        <h3 className="font-medium">Cardiology Follow-up</h3>
                                        <p className="text-sm text-gray-500">June 15, 2023 at 10:00 AM</p>
                                        <p className="text-sm text-gray-500">Dr. Robert Smith</p>
                                    </div>
                                    <Button>Reschedule</Button>
                                </div>
                                <div className="flex items-center space-x-4">
                                    <Avatar>
                                        <AvatarImage src="/placeholder.svg?height=40&width=40" alt="Dr. Sarah Johnson" />
                                        <AvatarFallback>SJ</AvatarFallback>
                                    </Avatar>
                                    <div className="flex-1">
                                        <h3 className="font-medium">Lab Results Review</h3>
                                        <p className="text-sm text-gray-500">June 22, 2023 at 2:00 PM</p>
                                        <p className="text-sm text-gray-500">Dr. Sarah Johnson</p>
                                    </div>
                                    <Button>Reschedule</Button>
                                </div>
                            </div>
                        </CardContent>
                    </Card>
                </TabsContent>
                <TabsContent value="medications">
                    <Card>
                        <CardHeader>
                            <CardTitle>Current Medications</CardTitle>
                        </CardHeader>
                        <CardContent>
                            <Table>
                                <TableHeader>
                                    <TableRow>
                                        <TableHead>Medication</TableHead>
                                        <TableHead>Dosage</TableHead>
                                        <TableHead>Frequency</TableHead>
                                        <TableHead>Prescribing Doctor</TableHead>
                                    </TableRow>
                                </TableHeader>
                                <TableBody>
                                    <TableRow>
                                        <TableCell>Lisinopril</TableCell>
                                        <TableCell>10mg</TableCell>
                                        <TableCell>Once daily</TableCell>
                                        <TableCell>Dr. Robert Smith</TableCell>
                                    </TableRow>
                                    <TableRow>
                                        <TableCell>Metformin</TableCell>
                                        <TableCell>500mg</TableCell>
                                        <TableCell>Twice daily</TableCell>
                                        <TableCell>Dr. Sarah Johnson</TableCell>
                                    </TableRow>
                                    <TableRow>
                                        <TableCell>Atorvastatin</TableCell>
                                        <TableCell>20mg</TableCell>
                                        <TableCell>Once daily at bedtime</TableCell>
                                        <TableCell>Dr. Robert Smith</TableCell>
                                    </TableRow>
                                </TableBody>
                            </Table>
                        </CardContent>
                    </Card>
                </TabsContent>
                <TabsContent value="notes">
                    <Card>
                        <CardHeader>
                            <CardTitle>Doctor's Notes</CardTitle>
                        </CardHeader>
                        <CardContent>
                            <div className="space-y-4">
                                <div className="bg-gray-50 p-4 rounded-lg">
                                    <div className="flex items-center space-x-2 mb-2">
                                        <Stethoscope className="text-blue-500" />
                                        <span className="font-medium">Dr. Sarah Johnson</span>
                                        <span className="text-sm text-gray-500">- May 15, 2023</span>
                                    </div>
                                    <p className="text-sm">
                                        Patient shows good progress in managing blood pressure. Recommend continuing current medication
                                        regimen and scheduling a follow-up in 3 months.
                                    </p>
                                </div>
                                <div className="bg-gray-50 p-4 rounded-lg">
                                    <div className="flex items-center space-x-2 mb-2">
                                        <Stethoscope className="text-blue-500" />
                                        <span className="font-medium">Dr. Robert Smith</span>
                                        <span className="text-sm text-gray-500">- March 22, 2023</span>
                                    </div>
                                    <p className="text-sm">
                                        Administered flu vaccine. Patient reported no immediate adverse reactions. Advised to monitor for
                                        any delayed symptoms and report if necessary.
                                    </p>
                                </div>
                            </div>
                        </CardContent>
                    </Card>
                </TabsContent>
            </Tabs>
        </>
    );
}

