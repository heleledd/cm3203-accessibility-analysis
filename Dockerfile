# Use an official lightweight Python image with Conda installed
FROM mambaorg/micromamba:1.5.1

# Set the working directory
WORKDIR /app

# Copy the environment definition (if you have an environment.yml)
# Or install dependencies directly. 
# It is safer to install osmnx and geopandas via conda than pip.
RUN micromamba install -y -n base -c conda-forge \
    python=3.9 \
    osmnx \
    geopandas \
    networkx \
    tqdm \
    scipy \
    && micromamba clean --all --yes

# Activate the environment
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Copy your source code
COPY src/ ./src/

# Create data directories (to be mounted as volumes)
RUN mkdir -p data/input_data data/output

# Set environment variables
ENV PYTHONPATH=/app/src
ENV NUM_WORKERS=4

# Run the application
CMD ["python", "src/main.py"]