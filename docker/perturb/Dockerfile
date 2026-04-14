FROM python:3.11-slim

RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    libhdf5-dev \
    pkg-config \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    scanpy==1.11.5 \
    anndata==0.10.7 \
    numpy==1.26.4 \
    pandas==2.2.2 \
    scipy==1.13.0 \
    scikit-learn==1.5.0 \
    leidenalg==0.10.2 \
    python-igraph==0.11.6

COPY bin/perturb_qc.py /usr/local/bin/perturb_qc.py
RUN chmod +x /usr/local/bin/perturb_qc.py

CMD ["python", "/usr/local/bin/perturb_qc.py"]