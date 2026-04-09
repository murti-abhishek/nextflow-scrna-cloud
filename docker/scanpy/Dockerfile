FROM python:3.11-slim

# System dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    libhdf5-dev \
    pkg-config \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies
RUN pip install --no-cache-dir \
    scanpy==1.10.1 \
    scrublet==0.2.3 \
    anndata==0.10.7 \
    numpy==1.26.4 \
    pandas==2.2.2 \
    scipy==1.13.0 \
    leidenalg==0.10.2 \
    python-igraph==0.11.6

# Copy QC script
COPY bin/scanpy_qc.py /usr/local/bin/scanpy_qc.py
RUN chmod +x /usr/local/bin/scanpy_qc.py

CMD ["python", "/usr/local/bin/scanpy_qc.py"]