"use client";

import { Button } from "@v1/ui/button";
import { useChat } from "ai/react";
import { revalidateTag } from "next/cache";

export function CreateCasesChat() {
  const { messages, input, handleInputChange, handleSubmit } = useChat({
    api: "/api/chat/cases",
    maxToolRoundtrips: 2,
  });
  return (
    <div className="flex flex-col w-full max-w-md py-24 mx-auto stretch relative">
      <div className="space-y-4 flex-grow">
        {messages.map((m) => (
          <div key={m.id} className="whitespace-pre-wrap">
            <div>
              <div className="font-bold">{m.role}</div>
              <p>
                {m.content.length > 0 ? (
                  m.content
                ) : (
                  <span className="italic font-light">
                    {`calling tool: ${m?.toolInvocations?.[0]?.toolName}`}
                  </span>
                )}
              </p>
            </div>
          </div>
        ))}
      </div>
      <form onSubmit={handleSubmit} className="absolute bottom-0 w-full">
        <input
          className="w-full p-2 mb-8 border border-gray-300 rounded shadow-xl"
          value={input}
          placeholder="Say something..."
          onChange={handleInputChange}
        />
      </form>
    </div>
  );
}
