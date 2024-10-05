// import { Header } from "@/components/header";
import { AppSidebar } from "@/components/app-sidebar";
import { getUser } from "@v1/supabase/queries";
import { SidebarLayout, SidebarTrigger } from "@v1/ui/sidebar";

export default async function Layout({
  children,
}: {
  children: React.ReactNode;
}) {
  const { cookies } = await import("next/headers");
  const { data: user } = await getUser();

  if (!user) {
    return <>No user found</>;
  }

  return (
    <SidebarLayout
      defaultOpen={cookies().get("sidebar:state")?.value === "true"}
    >
      <AppSidebar user={user} />
      <main className="flex flex-1 flex-col p-2 transition-all duration-300 ease-in-out">
        <div className="h-full p-2">
          <SidebarTrigger />
          {children}
        </div>
      </main>
    </SidebarLayout>
  );
}
