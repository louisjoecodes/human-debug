import { Database } from "lucide-react";

import { Card, CardContent } from "./card";
import { Progress } from "./progress";

export function StorageCard() {
  return (
    <Card className="rounded-md text-xs shadow-sm">
      <CardContent className="flex items-start gap-2.5 p-2.5">
        <div className="flex h-9 w-9 shrink-0 items-center justify-center rounded-md bg-accent text-accent-foreground">
          <Database className="h-5 w-5 text-muted-foreground" />
        </div>
        <div className="grid flex-1 gap-1">
          <p className="font-medium">Running out of AI credits?</p>
          <p className="text-muted-foreground">79.2K tokens / 11M used</p>
          <Progress
            value={7.2}
            className="mt-1"
            aria-label="79.2 GB / 100 GB used"
          />
        </div>
      </CardContent>
    </Card>
  );
}
