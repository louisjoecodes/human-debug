"use client";

import { findRelavantKnowledge } from "@v1/supabase/lib/ai/embedding";
import { Button } from "@v1/ui/button";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@v1/ui/select";
// import { useChat } from "ai/react";
// import { useState } from "react";

const systemMessages = {
  vet: "You are a friendly and professional veterinary receptionist. Use the knowledge base to answer questions about appointments, services, and pet care. If a caller wants to book an appointment, use the booking tool.",
  // Add more predefined system messages here
};

export default function PhoneBotSetup() {
  //   const [selectedRole, setSelectedRole] = useState("vet");

  //   const { messages, append, reload, stop, isLoading, input, setInput } =
  //     useChat({
  //       api: "/api/chat",
  //       initialMessages: [
  //         { role: "system", content: systemMessages[selectedRole] },
  //       ],
  //       body: {
  //         findRelavantKnowledge,
  //       },
  //     });

  //   const handleRoleChange = (value: string) => {
  //     setSelectedRole(value);
  //     reload([{ role: "system", content: systemMessages[value] }]);
  //   };

  return (
    <div className="container mx-auto p-4">
      <h1 className="text-2xl font-bold mb-4">Phone Bot Setup</h1>

      <div className="mb-4">
        <label className="block text-sm font-medium mb-2">
          Select Bot Role
        </label>
        <Select>
          <SelectTrigger className="w-full">
            <SelectValue placeholder="Select a role" />
          </SelectTrigger>
          <SelectContent>
            <SelectItem value="vet">Veterinary Receptionist</SelectItem>
            {/* Add more roles here */}
          </SelectContent>
        </Select>
      </div>

      <div className="mb-4">
        <h2 className="text-lg font-semibold mb-2">Chat Preview</h2>
        <div className="border p-4 h-64 overflow-y-auto">
          {/* {messages.map((m) => (
            <div key={m.id} className="mb-2">
              <strong>{m.role}:</strong> {m.content}
            </div>
          ))} */}
        </div>
      </div>

      <div className="mb-4">
        {/*
        <input
          className="w-full p-2 border rounded"
          value={input}
          onChange={(e) => setInput(e.target.value)}
          placeholder="Type a message to test the bot..."
        />
        <Button
          onClick={() => append({ role: "user", content: input })}
          disabled={isLoading}
          className="mt-2"
        >
          Send
        </Button> */}
      </div>

      <div className="mb-4">
        <h2 className="text-lg font-semibold mb-2">Twilio Integration</h2>
        <p>Twilio phone number: [Your Twilio number here]</p>
      </div>

      <div className="mb-4">
        <h2 className="text-lg font-semibold mb-2">ElevenLabs Integration</h2>
      </div>

      <Button className="mt-4">Save Configuration</Button>
    </div>
  );
}
