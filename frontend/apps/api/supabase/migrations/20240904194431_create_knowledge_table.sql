-- create knowledge table
create table knowledge (
  id uuid primary key default gen_random_uuid(),
  user_id uuid not null,
  content text not null,
  created_at timestamptz not null default now(),
  updated_at timestamptz not null default now()
);

-- add foreign key constraint
alter table
  knowledge
add
  constraint fk_knowledge_user foreign key (user_id) references public.users(id) on
delete
  cascade;

-- create index for faster queries
create index idx_knowledge_user_id on knowledge(user_id);

-- add rls policies
alter table
  knowledge enable row level security;

-- policy to allow read access for all authenticated users
create policy "allow read access for all authenticated users" on knowledge for
select
  to authenticated
  using (true);

-- policy to allow users to insert their own knowledge
create policy "allow insert for authenticated users" on knowledge for
insert
  with check (auth.uid() = user_id);

-- policy to allow users to update their own knowledge
create policy "allow update for knowledge owners" on knowledge for
update
  using (auth.uid() = user_id);

-- policy to allow users to delete their own knowledge
create policy "allow delete for knowledge owners" on knowledge for
delete
  using (auth.uid() = user_id);

-- trigger to call the update_updated_at function
create trigger update_knowledge_updated_at before
update
  on knowledge for each row execute function update_updated_at();