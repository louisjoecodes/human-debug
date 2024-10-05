"use client";

import { Button } from "@v1/ui/button";
import { RefreshCwIcon } from "lucide-react";
import { toast } from "sonner";

import { refreshCasesAction } from "@/actions/cases/refresh-cases-action";

export function CasesRefreshButton() {
  return (
    <Button
      variant="outline"
      size="icon"
      onClick={async () => {
        await refreshCasesAction({});
        toast.success("Cases refreshed ");
      }}
    >
      <RefreshCwIcon className="h-3 w-3" />
    </Button>
  );
}
