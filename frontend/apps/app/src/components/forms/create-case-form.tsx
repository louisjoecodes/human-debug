"use client";

import { createCaseAction } from "@/actions/cases/create-case-action";
import { createCaseSchema } from "@/actions/cases/schema";
import { zodResolver } from "@hookform/resolvers/zod";
import { Button } from "@v1/ui/button";
import {
  Form,
  FormControl,
  FormField,
  FormItem,
  FormMessage,
} from "@v1/ui/form";
import { Textarea } from "@v1/ui/textarea";
import { Input } from "@v1/ui/input";
import { useAction } from "next-safe-action/hooks";
import { useState } from "react";
import { useForm } from "react-hook-form";
import { toast } from "sonner";
import type { z } from "zod";
import { DatePickerDemo } from "@v1/ui/date-picker";
import { format } from "date-fns";

export function CreateCaseForm() {
  const form = useForm<z.infer<typeof createCaseSchema>>({
    resolver: zodResolver(createCaseSchema),
    defaultValues: {
      first_name: "",
      last_name: "",
      date_of_birth: "",
    },
  });
  const createCase = useAction(createCaseAction, {
    onExecute: () => { },
    onSuccess: () => {
      toast.success("Case created 🧠");
      form.reset({ first_name: "", last_name: "", date_of_birth: "" });
    },
    onError: () => { },
  });

  return (
    <Form {...form}>
      <form
        onSubmit={form.handleSubmit(createCase.execute)}
        className="space-y-4"
      >
        <FormField
          control={form.control}
          name="first_name"
          render={({ field }) => (
            <FormItem>
              <FormControl>
                <Input
                  placeholder="First name"
                  {...field}

                />
              </FormControl>
              <FormMessage />
            </FormItem>
          )}
        />
        <FormField
          control={form.control}
          name="last_name"
          render={({ field }) => (
            <FormItem>
              <FormControl>
                <Input
                  placeholder="Last name"
                  {...field}

                />
              </FormControl>
              <FormMessage />
            </FormItem>
          )}
        />
        <FormField
          control={form.control}
          name="date_of_birth"
          render={({ field }) => (
            <FormItem>
              <FormControl>
                <DatePickerDemo
                  date={field.value ? new Date(field.value) : undefined}
                  onSelect={(date) => field.onChange(date ? format(date, "yyyy-MM-dd") : "")}
                />
              </FormControl>
              <FormMessage />
            </FormItem>
          )}
        />
        <Button
          type="submit"
          disabled={createCase.status === "executing"}
          className="w-full"
        >
          {createCase.status === "executing"
            ? "Creating..."
            : "Create Knowledge"}
        </Button>
      </form>
    </Form>
  );
}
