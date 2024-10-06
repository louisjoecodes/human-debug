"use client";

import { useChat } from "ai/react";

export function AnalysisChat({
  patientData,
  phenotypes,
}: { patientData: any; phenotypes: any }) {
  const { messages, input, handleInputChange, handleSubmit } = useChat({
    api: "/api/chat/analysis?case_id=",
    maxToolRoundtrips: 2,
    initialMessages: [
      {
        id: "system-patient-data",
        role: "system",
        content: `content: ${JSON.stringify(patientData)}`,
      },
      {
        id: "system-phenotypes",
        role: "system",
        content: `content: ${JSON.stringify(phenotypes)}`,
      },
    ],
  });

  return (
    <div className="flex flex-col h-full">
      <div className="flex-grow overflow-y-auto" style={{ height: "calc(100% - 70px)" }}>
        <div className="space-y-4 p-4">
          {messages.map(
            (m) =>
              !m.content.startsWith("content:") && (
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
              )
          )}
        </div>
      </div>
      <form onSubmit={handleSubmit} className="p-4 border-t">
        <input
          className="w-full p-2 border border-gray-300 rounded shadow-xl"
          value={input}
          placeholder="Ask something..."
          onChange={handleInputChange}
        />
      </form>
    </div>
  );
}
