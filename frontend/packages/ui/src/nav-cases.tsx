import { type LucideIcon, MoreHorizontal, PlusSquare } from "lucide-react";
import Link from "next/link";

import { cn } from "@v1/ui";
import { Button } from "./button";
import { BriefcaseMedical } from "lucide-react";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuSeparator,
  DropdownMenuTrigger,
} from "./dropdown-menu";

export function NavCases({
  cases,
  className,
}: {
  cases: {
    id: string;
    first_name: string;
    last_name: string;
  }[];
} & React.ComponentProps<"ul">) {
  return (
    <ul className={cn("grid gap-0.5", className)}>
      {cases.map((caseItem) => (
        <li
          key={caseItem.first_name}
          className="has-[[data-state=open]]:bg-accent has-[[data-state=open]]:text-accent-foreground group relative rounded-md hover:bg-accent hover:text-accent-foreground"
        >
          <Link
            href={`/cases/${caseItem.id}`}
            className="flex h-7 items-center gap-2.5 overflow-hidden rounded-md px-1.5 text-xs outline-none ring-ring transition-all hover:bg-accent hover:text-accent-foreground focus-visible:ring-2"
          >
            <BriefcaseMedical className="h-4 w-4 shrink-0 translate-x-0.5 text-muted-foreground" />
            <div className="line-clamp-1 grow overflow-hidden pr-6 font-medium">
              {caseItem.first_name} {caseItem.last_name}
            </div>
          </Link>
          <DropdownMenu>
            <DropdownMenuTrigger asChild>
              <Button
                variant="ghost"
                className="peer absolute right-1 top-0.5 h-6 w-6 shrink-0 rounded-md bg-accent p-0 text-accent-foreground opacity-0 ring-ring transition-all focus-visible:ring-2 group-focus-within:opacity-100 group-hover:opacity-100 data-[state=open]:bg-accent data-[state=open]:opacity-100"
              >
                <MoreHorizontal className="h-4 w-4 text-muted-foreground" />
                <span className="sr-only">Toggle</span>
              </Button>
            </DropdownMenuTrigger>
            <DropdownMenuContent side="right" align="start" sideOffset={20}>
              <DropdownMenuItem>Share</DropdownMenuItem>
              <DropdownMenuItem>Rename</DropdownMenuItem>
              <DropdownMenuSeparator />
              <DropdownMenuItem>Archive</DropdownMenuItem>
            </DropdownMenuContent>
          </DropdownMenu>
        </li>
      ))}
      <li>
        <button className="flex h-7 w-full items-center gap-2.5 overflow-hidden rounded-md px-1.5 text-left text-xs ring-ring transition-all hover:bg-muted focus-visible:outline-none focus-visible:ring-2">
          <PlusSquare className="h-4 w-4 shrink-0 translate-x-0.5 text-muted-foreground" />
          <div className="line-clamp-1 overflow-hidden font-medium text-muted-foreground">
            New Case
          </div>
        </button>
      </li>
    </ul>
  );
}
