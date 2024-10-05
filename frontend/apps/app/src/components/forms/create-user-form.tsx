"use client";
import { createUserAction } from "@/actions/users/create-user-action";
import { createUserSchema } from "@/actions/users/schema";
import { zodResolver } from "@hookform/resolvers/zod";
import { Button } from "@v1/ui/button";
import { Input } from "@v1/ui/input";
import { useTransition } from "react";
import { useForm } from "react-hook-form";

type FormData = {
  full_name: string;
  email: string;
  role: string;
  phone: string;
};

export function CreateUserForm() {
  const [isPending, startTransition] = useTransition();
  const {
    register,
    handleSubmit,
    formState: { errors },
    reset,
  } = useForm<FormData>({
    resolver: zodResolver(createUserSchema),
  });

  const onSubmit = (data: FormData) => {
    startTransition(async () => {
      await createUserAction({ parsedInput: data });
      reset();
    });
  };

  return (
    <form onSubmit={handleSubmit(onSubmit)} className="space-y-4">
      <div>
        <Input {...register("full_name")} placeholder="Full Name" />
        {errors.full_name && (
          <p className="text-red-500 text-xs mt-1">
            {errors.full_name.message}
          </p>
        )}
      </div>
      <div>
        <Input {...register("email")} placeholder="Email" type="email" />
        {errors.email && (
          <p className="text-red-500 text-xs mt-1">{errors.email.message}</p>
        )}
      </div>
      <div>
        <Input {...register("role")} placeholder="Role" />
        {errors.role && (
          <p className="text-red-500 text-xs mt-1">{errors.role.message}</p>
        )}
      </div>
      <div>
        <Input {...register("phone")} placeholder="Phone" />
        {errors.phone && (
          <p className="text-red-500 text-xs mt-1">{errors.phone.message}</p>
        )}
      </div>
      <Button type="submit" disabled={isPending}>
        {isPending ? "Creating..." : "Create User"}
      </Button>
    </form>
  );
}
