"use client";

import {
  Atom,
  Frame,
  History,
  Hospital,
  LifeBuoy,
  Phone,
  PieChart,
  Send,
  BriefcaseMedical,
  SquareTerminal
} from "lucide-react";

import { NavUser } from "@/components/app-sidebar-nav-user";
import type { Database } from "@v1/supabase/types";
import { NavMain } from "@v1/ui/nav-main";
import { NavCases } from "@v1/ui/nav-cases";
import { NavSecondary } from "@v1/ui/nav-secondary";
import {
  Sidebar,
  SidebarContent,
  SidebarFooter,
  SidebarHeader,
  SidebarItem,
  SidebarLabel,
} from "@v1/ui/sidebar";
import { StorageCard } from "@v1/ui/storage-card";
import { TeamSwitcher } from "@v1/ui/team-switcher";
import { createClient } from "@v1/supabase/client";
import { useEffect, useState } from "react";
const data = {
  teams: [
    {
      name: "Genomics Lab 🧪",
      logo: Atom,
      plan: "Enterprise",
    }
  ],
  user: {
    name: "shadcn",
    email: "m@example.com",
    avatar: "/avatars/shadcn.jpg",
  },
  navMain: [
    {
      title: "Case Management",
      url: "#",
      icon: Hospital,
      isActive: true,
      items: [
        {
          title: "Cases",
          url: "/cases",
          icon: History,
          description: "View your cases",
        }
      ],
    },
  ],

  navSecondary: [
    {
      title: "Support",
      url: "#",
      icon: LifeBuoy,
    },
    {
      title: "Feedback",
      url: "#",
      icon: Send,
    },
  ],
  searchResults: [
    {
      title: "Routing Fundamentals",
      teaser:
        "The skeleton of every application is routing. This page will introduce you to the fundamental concepts of routing for the web and how to handle routing in Next.js.",
      url: "#",
    },
    {
      title: "Layouts and Templates",
      teaser:
        "The special files layout.js and template.js allow you to create UI that is shared between routes. This page will guide you through how and when to use these special files.",
      url: "#",
    },
    {
      title: "Data Fetching, Caching, and Revalidating",
      teaser:
        "Data fetching is a core part of any application. This page goes through how you can fetch, cache, and revalidate data in React and Next.js.",
      url: "#",
    },
    {
      title: "Server and Client Composition Patterns",
      teaser:
        "When building React applications, you will need to consider what parts of your application should be rendered on the server or the client. ",
      url: "#",
    },
    {
      title: "Server Actions and Mutations",
      teaser:
        "Server Actions are asynchronous functions that are executed on the server. They can be used in Server and Client Components to handle form submissions and data mutations in Next.js applications.",
      url: "#",
    },
  ],
};

export function AppSidebar({
  user,
}: { user: Database["public"]["Tables"]["users"]["Row"] }) {
  const client = createClient();

  const [cases, setCases] = useState([]);

  useEffect(() => {
    const fetchCases = async () => {
      const { data, error } = await client.from("cases").select("*");
      if (error) {
        console.error("Error fetching cases:", error);
      } else {
        setCases(data as any);
      }
    };

    fetchCases();
  }, []);

  console.log(cases);
  console.log("____")
  return (
    <Sidebar>
      <SidebarHeader>
        <TeamSwitcher teams={data.teams} />
      </SidebarHeader>
      <SidebarContent>
        <SidebarItem>
          <SidebarLabel>Platform</SidebarLabel>
          <NavMain items={data.navMain} searchResults={data.searchResults} />
        </SidebarItem>
        <SidebarItem>
          <SidebarLabel>Cases</SidebarLabel>
          <NavCases cases={cases} />
        </SidebarItem>
        <SidebarItem className="mt-auto">
          <SidebarLabel>Help</SidebarLabel>
          <NavSecondary items={data.navSecondary} />
        </SidebarItem>
        <SidebarItem>
          <StorageCard />
        </SidebarItem>
      </SidebarContent>
      <SidebarFooter>
        <NavUser user={user} />
      </SidebarFooter>
    </Sidebar>
  );
}
