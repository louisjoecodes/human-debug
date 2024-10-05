import { CreateCasesChat } from "@/components/create-cases-chat";
import { CreateCaseForm } from "@/components/forms/create-case-form";
import { CasesRefreshButton } from "@/components/cases/cases-refresh-button";

import { CasesLoading } from "@/components/cases/cases.loading";
import { CasesServer } from "@/components/cases/cases.server";
import {
    Card,
    CardContent,
    CardFooter,
    CardHeader,
    CardTitle,
} from "@v1/ui/card";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@v1/ui/tabs";
import { Suspense } from "react";

export const metadata = {
    title: "Cases",
};

export default function Page() {
    return (
        <>
            <div className="p-4 sm:px-6 sm:py-0 mb-3">
                <h1 className="text-lg font-semibold md:text-2xl">Cases</h1>
            </div>

            <div className="grid flex-1 items-start gap-4 p-4 sm:px-6 sm:py-0 md:gap-8 lg:grid-cols-3 xl:grid-cols-3">
                <div className="grid auto-rows-max items-start gap-4 md:gap-8 lg:col-span-2">

                    <Suspense fallback={<CasesLoading />}>
                        <CasesServer />
                    </Suspense>
                </div>
                <div>
                    <Card className="overflow-hidden">
                        <CardHeader className="flex flex-row items-start bg-muted/50">
                            <div className="grid gap-0.5">
                                <CardTitle className="group flex items-center gap-2 text-lg">
                                    Create Case
                                </CardTitle>
                            </div>
                        </CardHeader>
                        <CardContent className="p-6 text-sm">
                            <div className="grid gap-3">
                                <Tabs defaultValue="manual">
                                    <TabsList>
                                        <TabsTrigger value="chat">Chat</TabsTrigger>
                                        <TabsTrigger value="manual">Manual</TabsTrigger>
                                    </TabsList>
                                    <TabsContent value="chat">
                                        <CreateCasesChat />
                                    </TabsContent>
                                    <TabsContent value="manual">
                                        <CreateCaseForm />
                                    </TabsContent>
                                </Tabs>
                            </div>
                        </CardContent>
                        <CardFooter className="flex flex-row items-center border-t bg-muted/50 px-6 py-3">
                            <div className="text-xs text-muted-foreground">Experimental</div>
                        </CardFooter>
                    </Card>
                </div>
            </div>
        </>
    );
}
