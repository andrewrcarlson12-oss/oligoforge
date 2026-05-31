# OligoForge as an always-updated web app.
# Any host that builds Dockerfiles (Render, Railway, Fly.io, a VPS) can run this,
# and redeploy automatically whenever you push to GitHub.
FROM python:3.12-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
# NCBI requires an email; set real values as host environment variables.
ENV OLIGOFORGE_EMAIL=arcarl27@colby.edu
# Optional NCBI API key (10 req/s vs 3). Pass at runtime, never bake a secret into the image:
#   docker run -e OLIGOFORGE_NCBI_KEY=your_key ...
# Hosts inject $PORT; fall back to 8111 locally.
CMD ["sh", "-c", "uvicorn app:app --host 0.0.0.0 --port ${PORT:-8111} --limit-concurrency 24 --timeout-keep-alive 30"]
