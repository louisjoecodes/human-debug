CREATE EXTENSION IF NOT EXISTS vector;

-- Create knowledge_content_embeddings table
CREATE TABLE knowledge_content_embeddings (
  id uuid primary key default gen_random_uuid(),
  knowledge_id UUID REFERENCES knowledge(id) ON DELETE CASCADE,
  content TEXT NOT NULL,
  embedding VECTOR(1536) NOT NULL
);

-- Create index using HNSW for vector similarity search
CREATE INDEX knowledge_content_embeddingIndex ON knowledge_content_embeddings USING hnsw (embedding vector_cosine_ops);

-- Enable row level security
ALTER TABLE knowledge_content_embeddings ENABLE ROW LEVEL SECURITY;

-- Policy to allow read access for all authenticated users
CREATE POLICY "Allow authenticated users to read knowledge content embeddings" ON knowledge_content_embeddings
  FOR SELECT TO authenticated USING (true);

-- Policy to allow insert for authenticated users who own the related knowledge
CREATE POLICY "Allow knowledge owners to insert content embeddings" ON knowledge_content_embeddings
  FOR INSERT WITH CHECK (
    EXISTS (
      SELECT 1 FROM knowledge
      WHERE knowledge.id = knowledge_content_embeddings.knowledge_id
        AND knowledge.user_id = auth.uid()
    )
  );

-- Policy to allow update for knowledge owners
CREATE POLICY "Allow knowledge owners to update content embeddings" ON knowledge_content_embeddings
  FOR UPDATE USING (
    EXISTS (
      SELECT 1 FROM knowledge
      WHERE knowledge.id = knowledge_content_embeddings.knowledge_id
        AND knowledge.user_id = auth.uid()
    )
  );

-- Policy to allow delete for knowledge owners
CREATE POLICY "Allow knowledge owners to delete content embeddings" ON knowledge_content_embeddings
  FOR DELETE USING (
    EXISTS (
      SELECT 1 FROM knowledge
      WHERE knowledge.id = knowledge_content_embeddings.knowledge_id
        AND knowledge.user_id = auth.uid()
    )
  );

-- Create trigger to update updated_at in knowledge table
CREATE OR REPLACE FUNCTION update_knowledge_updated_at()
RETURNS TRIGGER AS $$
BEGIN
  UPDATE knowledge
  SET updated_at = NOW()
  WHERE id = NEW.knowledge_id;
  RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER update_knowledge_updated_at_on_embedding_change
AFTER INSERT OR UPDATE OR DELETE ON knowledge_content_embeddings
FOR EACH ROW EXECUTE FUNCTION update_knowledge_updated_at();

-- Create cosine similarity search function
CREATE OR REPLACE FUNCTION cosine_similarity_search_knowledge(query_embedding VECTOR(1536))
RETURNS TABLE (content TEXT, similarity FLOAT) AS $$
BEGIN
  RETURN QUERY
  SELECT e.content, 1 - (e.embedding <=> query_embedding) AS similarity
  FROM knowledge_content_embeddings e
  WHERE 1 - (e.embedding <=> query_embedding) > 0.5 -- Only return similar items
  ORDER BY similarity DESC
  LIMIT 4;
END;
$$ LANGUAGE plpgsql;