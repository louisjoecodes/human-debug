import { CreateUserForm } from "@/components/forms/create-user-form";
import { UserList } from "@/components/users/user-list";
import {
  Card,
  CardContent,
  CardFooter,
  CardHeader,
  CardTitle,
} from "@v1/ui/card";
import { Suspense } from "react";

export const metadata = {
  title: "Users",
};

export default function Page() {
  return (
    <>
      <div className="p-4 sm:px-6 sm:py-0 mb-3">
        <h1 className="text-lg font-semibold md:text-2xl">Users</h1>
      </div>

      <div className="grid flex-1 items-start gap-4 p-4 sm:px-6 sm:py-0 md:gap-8 lg:grid-cols-3 xl:grid-cols-3">
        <div className="grid auto-rows-max items-start gap-4 md:gap-8 lg:col-span-2">
          <div className="flex justify-end">{/* <UserRefreshButton /> */}</div>
          <Suspense fallback={<div>Loading...</div>}>
            <UserList />
          </Suspense>
        </div>
        <div>
          <Card className="overflow-hidden">
            <CardHeader className="flex flex-row items-start bg-muted/50">
              <div className="grid gap-0.5">
                <CardTitle className="group flex items-center gap-2 text-lg">
                  Add User
                </CardTitle>
              </div>
            </CardHeader>
            <CardContent className="p-6 text-sm">
              <CreateUserForm />
            </CardContent>
            <CardFooter className="flex flex-row items-center border-t bg-muted/50 px-6 py-3">
              <div className="text-xs text-muted-foreground">Manage users</div>
            </CardFooter>
          </Card>
        </div>
      </div>
    </>
  );
}
