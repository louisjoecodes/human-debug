"use client";

import { deleteKnowledgeAction } from "@/actions/knowledge/delete-knowledge-action";
import { Button } from "@v1/ui/button";
import { TrashIcon } from "lucide-react";
import { toast } from "sonner";

export function KnowledgeDeleteButton({ id }: { id: string }) {
  return (
    <Button
      variant="outline"
      size="icon"
      onClick={async () => {
        await deleteKnowledgeAction({ id });
        toast.success("Knowledge deleted 🗑️");
      }}
    >
      <TrashIcon className="h-3 w-3" />
    </Button>
  );
}
