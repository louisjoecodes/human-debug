-- Create members_role_description_embeddings table
CREATE TABLE members_role_description_embeddings (
  id uuid primary key default gen_random_uuid(),
  member_id UUID REFERENCES members(id) ON DELETE CASCADE,
  role_description TEXT NOT NULL,
  embedding VECTOR(1536) NOT NULL
);

-- Create index using HNSW for vector similarity search
CREATE INDEX member_role_description_embeddingIndex ON members_role_description_embeddings USING hnsw (embedding vector_cosine_ops);

-- Enable row level security
ALTER TABLE members_role_description_embeddings ENABLE ROW LEVEL SECURITY;

-- Policy to allow read access for all authenticated users
CREATE POLICY "Allow authenticated users to read member role_description embeddings" ON members_role_description_embeddings
  FOR SELECT TO authenticated USING (true);

-- Policy to allow insert for authenticated users who own the related member
CREATE POLICY "Allow member owners to insert role_description embeddings" ON members_role_description_embeddings
  FOR INSERT WITH CHECK (
    EXISTS (
      SELECT 1 FROM members
      WHERE members.id = members_role_description_embeddings.member_id
        AND members.user_id = auth.uid()
    )
  );

-- Policy to allow update for member owners
CREATE POLICY "Allow member owners to update content embeddings" ON members_role_description_embeddings
  FOR UPDATE USING (
    EXISTS (
      SELECT 1 FROM members
      WHERE members.id = members_role_description_embeddings.member_id
        AND members.user_id = auth.uid()
    )
  );

-- Policy to allow delete for member owners
CREATE POLICY "Allow knowledge owners to delete member embeddings" ON members_role_description_embeddings
  FOR DELETE USING (
    EXISTS (
      SELECT 1 FROM members
      WHERE members.id = members_role_description_embeddings.member_id
        AND members.user_id = auth.uid()
    )
  );

-- Create trigger to update updated_at in members table
CREATE OR REPLACE FUNCTION update_members_updated_at()
RETURNS TRIGGER AS $$
BEGIN
  UPDATE members
  SET updated_at = NOW()
  WHERE id = NEW.member_id;
  RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER update_member_updated_at_on_embedding_change
AFTER INSERT OR UPDATE OR DELETE ON members_role_description_embeddings
FOR EACH ROW EXECUTE FUNCTION update_members_updated_at();

-- Create cosine similarity search function
CREATE OR REPLACE FUNCTION cosine_similarity_search_members(query_embedding VECTOR(1536))
RETURNS TABLE (role_description TEXT, similarity FLOAT) AS $$
BEGIN
  RETURN QUERY
  SELECT e.content, 1 - (e.embedding <=> query_embedding) AS similarity
  FROM members_role_description_embeddings e
  WHERE 1 - (e.embedding <=> query_embedding) > 0.5 -- Only return similar items
  ORDER BY similarity DESC
  LIMIT 4;
END;
$$ LANGUAGE plpgsql;