FROM python:3.11-slim

# ── system dependencies ──────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    minimap2 \
    ncbi-blast+ \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /viralift

# ── python dependencies ──────────────────────────────────────────
COPY ui/requirements.txt /viralift/ui/requirements.txt
RUN pip install --no-cache-dir -r ui/requirements.txt

# ── application code ─────────────────────────────────────────────
COPY app  /viralift/app
COPY ui   /viralift/ui

# ── runtime ──────────────────────────────────────────────────────
ENV PYTHONPATH=/viralift
EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health || exit 1

CMD ["streamlit", "run", "ui/streamlit_app.py", \
     "--server.port=8501", \
     "--server.address=0.0.0.0", \
     "--server.headless=true"]