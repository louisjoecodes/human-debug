-- create case table
create table cases (
  id uuid primary key default gen_random_uuid(),
  user_id uuid not null,
  first_name text not null,
  last_name text not null,
  date_of_birth date not null,
  created_at timestamptz not null default now(),
  updated_at timestamptz not null default now()
);

-- add foreign key constraint
alter table
  cases
add
  constraint fk_cases_user foreign key (user_id) references public.users(id) on
delete
  cascade;

-- create index for faster queries
create index idx_cases_user_id on cases(user_id);

-- add rls policies
alter table
  cases enable row level security;

-- policy to allow read access for all authenticated users
create policy "allow read access for all authenticated users" on cases for
select
  to authenticated
  using (true);

-- policy to allow users to insert their own cases
create policy "allow insert for authenticated users" on cases for
insert
  with check (auth.uid() = user_id);

-- policy to allow users to update their own cases
create policy "allow update for cases owners" on cases for
update
  using (auth.uid() = user_id);

-- policy to allow users to delete their own cases
create policy "allow delete for cases owners" on cases for
delete
  using (auth.uid() = user_id);

-- trigger to call the update_updated_at function
create trigger update_cases_updated_at before
update
  on cases for each row execute function update_updated_at();