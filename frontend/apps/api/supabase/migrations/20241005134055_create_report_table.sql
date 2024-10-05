-- Step 1: Create the reports table
CREATE TABLE reports (
  id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
  case_id uuid NOT NULL,
  created_at timestamptz NOT NULL DEFAULT NOW(),
  updated_at timestamptz NOT NULL DEFAULT NOW()
);

-- Add foreign key constraint
ALTER TABLE reports
ADD CONSTRAINT fk_reports_case FOREIGN KEY (case_id)
REFERENCES public.cases(id) ON DELETE CASCADE;

-- Create index on case_id
CREATE INDEX idx_reports_cases_id ON reports(case_id);

-- Create trigger for updated_at
CREATE TRIGGER update_reports_updated_at
BEFORE UPDATE ON reports
FOR EACH ROW EXECUTE FUNCTION update_updated_at();

-- Step 2: Enable Row-Level Security (Optional)
ALTER TABLE reports ENABLE ROW LEVEL SECURITY;

-- Allow full access to reports
CREATE POLICY "allow full access to reports" ON reports
USING (true)
WITH CHECK (true);

-- Step 3: Create a public storage bucket named 'reports'
INSERT INTO storage.buckets
  (id, name)
VALUES
  ('reports', 'reports');

-- Step 4: Configure storage policies
-- Public read access to the 'reports' bucket
CREATE POLICY "Public read access to reports bucket" 
ON storage.objects FOR SELECT
USING (bucket_id = 'reports');

-- Authenticated write access to the 'reports' bucket
CREATE POLICY "Authenticated write access to reports bucket" 
ON storage.objects FOR INSERT
WITH CHECK (
  bucket_id = 'reports' 
  AND auth.role() = 'authenticated'
);

-- Authenticated update access to the 'reports' bucket
CREATE POLICY "Authenticated update access to reports bucket" 
ON storage.objects FOR UPDATE
USING (
  bucket_id = 'reports' 
  AND auth.role() = 'authenticated'
)
WITH CHECK (bucket_id = 'reports');

-- Authenticated delete access to the 'reports' bucket
CREATE POLICY "Authenticated delete access to reports bucket" 
ON storage.objects FOR DELETE
USING (
  bucket_id = 'reports' 
  AND auth.role() = 'authenticated'
);
