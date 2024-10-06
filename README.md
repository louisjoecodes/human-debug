
![Frame 27](https://github.com/user-attachments/assets/8a118387-3c94-4d28-a97e-acc088a6ad75)



# Human Debug 🧬

Human Debug is an advanced genomic analysis platform designed to transform how clinicians interpret genetic data and manage patient cases, facilitating precision medicine.

## 🚀 Features

- **Comprehensive Genomic Analysis**: Advanced algorithms for in-depth genomic data interpretation.
- **Patient Case Management**: Streamlined management and tracking of patient cases for improved clinical workflow.
- **Interactive Genomic Visualizations**: Heatmaps and other visual tools for enhanced data interpretation.
- **AI-Powered Genomic Insights**: Leverage artificial intelligence for advanced genetic variant analysis and interpretation.

## 🛠️ Tech Stack

- **Frontend**: Next.js, React, TypeScript
- **Backend**: Python (FastAPI)
- **Database**: Supabase
- **UI Components**: Shadcn UI
- **Styling**: Tailwind CSS
- **State Management**: React Hooks
- **Authentication**: Supabase Auth
- **File Handling**: React Dropzone
- **Visualization**: React Heatmap Grid, AG Charts
- **AI Integration**: OpenAI

## 🚀 Getting Started

### Prerequisites

- Node.js
- Python 3.8+
- Docker
- Supabase account

### Frontend Setup

1. Navigate to the frontend directory:

   ```bash
   cd frontend
   ```

2. Install dependencies:

   ```bash
   bun install
   ```

3. Set up environment variables:

   ```bash
   cp apps/app/.env.example apps/app/.env
   ```

   Edit the `.env` file with your Supabase and other API keys.

4. Run the development server:
   ```bash
   bun dev
   ```

### Backend Setup

1. Navigate to the backend directory:

   ```bash
   cd backend
   ```

2. Create a virtual environment and activate it:

   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
   ```

3. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

4. Run the FastAPI server:
   ```bash
   uvicorn src.main:app --reload
   ```

## 📚 Documentation

For more detailed information on how to use and contribute to Human Debug, please refer to our [documentation](link-to-your-docs).

## 🙏 Acknowledgements

- [Supabase](https://supabase.io/)
- [Next.js](https://nextjs.org/)
- [FastAPI](https://fastapi.tiangolo.com/)
- [Shadcn UI](https://ui.shadcn.com/)

---
