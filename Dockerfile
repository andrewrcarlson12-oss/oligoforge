# OligoForge web service. Render/Railway/Fly.io can rebuild this image on each push.
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .

# Do not run a public scientific service as root.
RUN useradd --create-home --uid 10001 oligoforge \
    && chown -R oligoforge:oligoforge /app
USER oligoforge

# Set OLIGOFORGE_EMAIL and optional OLIGOFORGE_NCBI_KEY at deployment time.
# Hosts inject $PORT; fall back to 8111 locally.
CMD ["sh", "-c", "uvicorn app:app --host 0.0.0.0 --port ${PORT:-8111} --limit-concurrency 24 --timeout-keep-alive 30 --no-proxy-headers"]
