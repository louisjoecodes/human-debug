"use client";

import { createKnowledgeAction } from "@/actions/knowledge/create-knowledge-action";
import { createKnowledgeSchema } from "@/actions/knowledge/schema";
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
import { useAction } from "next-safe-action/hooks";
import { useState } from "react";
import { useForm } from "react-hook-form";
import { toast } from "sonner";
import type { z } from "zod";

export function CreateKnowledgeForm() {
  const form = useForm<z.infer<typeof createKnowledgeSchema>>({
    resolver: zodResolver(createKnowledgeSchema),
    defaultValues: {
      content: "",
    },
  });
  const createKnowledge = useAction(createKnowledgeAction, {
    onExecute: () => {},
    onSuccess: () => {
      toast.success("Knowledge created 🧠");
      form.reset({ content: "" });
    },
    onError: () => {},
  });

  return (
    <Form {...form}>
      <form
        onSubmit={form.handleSubmit(createKnowledge.execute)}
        className="space-y-4"
      >
        <FormField
          control={form.control}
          name="content"
          render={({ field }) => (
            <FormItem>
              <FormControl>
                <Textarea
                  placeholder="Manually add your knowledge here."
                  {...field}
                  className="min-h-[120px] text-base"
                />
              </FormControl>
              <FormMessage />
            </FormItem>
          )}
        />
        <Button
          type="submit"
          disabled={createKnowledge.status === "executing"}
          className="w-full"
        >
          {createKnowledge.status === "executing"
            ? "Creating..."
            : "Create Knowledge"}
        </Button>
      </form>
    </Form>
  );
}
