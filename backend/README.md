### Human Debug Backend (API)

#### Setup

1. Install uv if you haven't already:
   ```
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

2. Install required system dependencies (for macOS users):
   ```
   brew install freebayes
   brew install bwa
   brew install samtools
   brew install bcftools
   ```

3. Create a virtual environment and install dependencies:
   ```
   cd backend
   uv venv
   source .venv/bin/activate  # On Unix-like systems
   ```

4. Start the development server:
   ```
   uv run src.main:start
   ```

The API should now be running at `http://localhost:8000`.

#### Development

- To add new dependencies:
  ```
  uv add <package_name>
  ```

- To update dependencies:
  ```
  uv pip compile pyproject.toml -o requirements.txt
  uv pip sync requirements.txt
  ```

- To run tests (assuming pytest is set up):
  ```
  uv run pytest
  ```
