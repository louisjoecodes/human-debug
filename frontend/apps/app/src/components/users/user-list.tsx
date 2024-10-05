import { getUsers } from "@/actions/users/get-users-action";

export async function UserList() {
  const { data: users, error } = await getUsers();

  if (error) {
    return <div>Error loading users: {error.message}</div>;
  }

  return (
    <div className="grid gap-4">
      {users?.map((user) => (
        <div key={user.id} className="bg-card p-4 rounded-lg shadow">
          <h3 className="font-semibold">{user.full_name}</h3>
          <p className="text-sm text-muted-foreground">
            {user.role || "No role specified"}
          </p>
          <p className="text-sm">{user.phone || "No phone number"}</p>
        </div>
      ))}
    </div>
  );
}
