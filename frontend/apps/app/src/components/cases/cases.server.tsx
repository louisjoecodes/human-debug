import { CaseDeleteButton } from "@/components/cases/cases-delete-button";
import Link from "next/link";
import { getCases } from "@v1/supabase/queries";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@v1/ui/table";
import { ScrollArea } from "@v1/ui/scroll-area";
import { formatDistanceToNow } from "date-fns";
import { Button } from "@v1/ui/button";
import { EyeIcon, PlusCircle } from "lucide-react";
import { Badge } from "@v1/ui/badge";
import { Card, CardContent, CardDescription, CardFooter, CardHeader, CardTitle } from "@v1/ui/card";
import { CasesRefreshButton } from "@/components/cases/cases-refresh-button";

export async function CasesServer() {
  const { data } = await getCases();

  return (
    <Card>
      <CardHeader>

        <div className="flex justify-between items-center">
          <CardTitle>Patient Cases</CardTitle>
          <CasesRefreshButton />
        </div>

        <CardDescription>
          Manage patient cases and view case details.
        </CardDescription>

      </CardHeader>
      <CardContent>

        <div className="flex justify-between items-center mb-4">
          <div className="text-sm text-muted-foreground">
            Showing <strong>{data?.length || 0}</strong> cases
          </div>

        </div>
        <ScrollArea className="h-[calc(100vh-300px)]">
          <Table>
            <TableHeader>
              <TableRow>
                <TableHead>First Name</TableHead>
                <TableHead>Last Name</TableHead>
                <TableHead>Birthdate</TableHead>
                <TableHead>Created</TableHead>
                <TableHead className="text-right">Actions</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {data
                ?.sort(
                  (a, b) =>
                    new Date(b.created_at).getTime() -
                    new Date(a.created_at).getTime(),
                )
                .map((caseItem) => (
                  <TableRow key={caseItem.id}>
                    <TableCell>{caseItem.first_name}</TableCell>
                    <TableCell>{caseItem.last_name}</TableCell>
                    <TableCell>{caseItem.date_of_birth}</TableCell>
                    <TableCell>
                      {formatDistanceToNow(new Date(caseItem.created_at), {
                        addSuffix: true,
                      })}
                    </TableCell>
                    <TableCell className="text-right space-x-2">
                      <Link href={`/cases/${caseItem.id}`} passHref>
                        <Button variant="outline" size="sm">
                          <EyeIcon className="h-4 w-4 mr-2" />
                          View
                        </Button>
                      </Link>
                      <CaseDeleteButton id={caseItem.id} />
                    </TableCell>
                  </TableRow>
                ))}
            </TableBody>
          </Table>
        </ScrollArea>
      </CardContent>
    </Card>
  );
}
