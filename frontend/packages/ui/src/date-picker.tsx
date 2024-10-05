"use client";

import * as React from "react";
import { format } from "date-fns";
import { Calendar as CalendarIcon } from "lucide-react";

import { cn } from "@v1/ui";
import { Button } from "@v1/ui/button";
import { Calendar } from "@v1/ui/calendar";
import { Popover, PopoverContent, PopoverTrigger } from "@v1/ui/popover";

interface DatePickerDemoProps {
  date?: Date;
  onSelect: (date: Date | undefined) => void;
}

export function DatePickerDemo({ date, onSelect }: DatePickerDemoProps) {
  return (
    <Popover>
      <PopoverTrigger asChild>
        <Button
          variant={"outline"}
          className={cn(
            "w-full justify-start text-left font-normal",
            !date && "text-muted-foreground",
          )}
        >
          <CalendarIcon className="mr-2 h-4 w-4" />
          {date ? format(date, "PPP") : <span>Pick a date</span>}
        </Button>
      </PopoverTrigger>
      <PopoverContent className="w-auto p-0">
        <Calendar
          mode="single"
          selected={date}
          onSelect={onSelect}
          initialFocus
        />
      </PopoverContent>
    </Popover>
  );
}
