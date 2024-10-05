-- create members table
create table members (
  id uuid primary key default gen_random_uuid(),
  user_id uuid not null,
  full_name text not null,
  phone_number text not null,
  role_description text not null,
  created_at timestamptz not null default now(),
  updated_at timestamptz not null default now()
);

-- add foreign key constraint
alter table
  members
add
  constraint fk_members_user foreign key (user_id) references public.users(id) on
delete
  cascade;

-- create index for faster queries
create index idx_members_user_id on members(user_id);

-- add rls policies
alter table
  members enable row level security;

-- policy to allow read access for all authenticated users
create policy "allow read access for all authenticated users" on members for
select
  to authenticated
  using (true);

-- policy to allow users to insert their own members
create policy "allow insert for authenticated users" on members for
insert
  with check (auth.uid() = user_id);

-- policy to allow users to update their own members
create policy "allow update for member owners" on members for
update
  using (auth.uid() = user_id);

-- policy to allow users to delete their own members
create policy "allow delete for member owners" on members for
delete
  using (auth.uid() = user_id);

-- trigger to call the update_updated_at function
create trigger update_members_updated_at before
update
  on members for each row execute function update_updated_at();