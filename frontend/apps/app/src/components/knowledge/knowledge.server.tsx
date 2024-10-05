import { KnowledgeDeleteButton } from "@/components/knowledge/knowledge-delete-button";

import { getKnowledge } from "@v1/supabase/queries";
import { Card, CardDescription, CardHeader } from "@v1/ui/card";
import { ScrollArea } from "@v1/ui/scroll-area";
import { formatDistanceToNow } from "date-fns";
import { LightbulbIcon } from "lucide-react";

export async function KnowledegServer() {
  const { data } = await getKnowledge();
  return (
    <ScrollArea className="h-[calc(100vh)] pr-4">
      <div className="grid gap-6 sm:grid-cols-2 lg:grid-cols-3">
        {data
          ?.sort(
            (a, b) =>
              new Date(b.created_at).getTime() -
              new Date(a.created_at).getTime(),
          )
          .map((knowledge) => (
            <Card key={knowledge.id}>
              <CardHeader>
                <div className="flex flex-row justify-between items-center">
                  <LightbulbIcon className={"h-3 w-3"} />
                  <div className="flex items-center gap-2">
                    <div className="text-sm text-muted-foreground">
                      {formatDistanceToNow(new Date(knowledge.created_at), {
                        addSuffix: true,
                      })}
                    </div>
                    <KnowledgeDeleteButton id={knowledge.id} />
                  </div>
                </div>
                <CardDescription>{knowledge.content}</CardDescription>
              </CardHeader>
            </Card>
          ))}
      </div>
    </ScrollArea>
  );
}
