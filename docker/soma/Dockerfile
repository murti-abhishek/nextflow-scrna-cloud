FROM python:3.11-slim

RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    pkg-config \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    tiledbsoma==2.3.0 \
    anndata==0.11.2 \
    pandas==2.2.2 \
    numpy==1.26.4 \
    more-itertools

COPY bin/ingest_soma.py /usr/local/bin/ingest_soma.py
RUN chmod +x /usr/local/bin/ingest_soma.py

CMD ["python", "/usr/local/bin/ingest_soma.py"]