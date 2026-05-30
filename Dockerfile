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
# Hosts inject $PORT; fall back to 8111 locally.
CMD ["sh", "-c", "uvicorn app:app --host 0.0.0.0 --port ${PORT:-8111}"]
