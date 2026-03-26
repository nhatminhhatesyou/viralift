FROM python:3.11-slim

RUN apt-get update && apt-get install -y \
    minimap2 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY app/requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY app /app

ENTRYPOINT ["python", "-m", "src.main"]