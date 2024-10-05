"use client";

import { Button } from "@v1/ui/button";
import { RefreshCwIcon } from "lucide-react";
import { toast } from "sonner";

import { refreshKnowledgeAction } from "@/actions/knowledge/refresh-knowledge-action";

export function KnowledgeRefreshButton() {
  return (
    <Button
      variant="outline"
      size="icon"
      onClick={async () => {
        await refreshKnowledgeAction({});
        toast.success("Knowledge refreshed ");
      }}
    >
      <RefreshCwIcon className="h-3 w-3" />
    </Button>
  );
}
