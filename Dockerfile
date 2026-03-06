FROM python:3.10-slim

# Install system dependencies required for geospatial libraries (GDAL, etc.)
RUN apt-get update && apt-get install -y \
    gdal-bin \
    libgdal-dev \
    g++ \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy source code
COPY src/ ./src/

# Set python path
ENV PYTHONPATH=/app

CMD ["python", "src/grid_to_features.py"]