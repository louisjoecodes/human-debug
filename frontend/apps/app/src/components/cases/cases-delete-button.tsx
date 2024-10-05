"use client";

import { deleteCaseAction } from "@/actions/cases/delete-case-action";
import { Button } from "@v1/ui/button";
import { TrashIcon } from "lucide-react";
import { toast } from "sonner";

export function CaseDeleteButton({ id }: { id: string }) {
  return (
    <Button
      variant="outline"
      size="icon"
      onClick={async () => {
        await deleteCaseAction({ id });
        toast.success("Case deleted 🗑️");
      }}
    >
      <TrashIcon className="h-3 w-3" />
    </Button>
  );
}
