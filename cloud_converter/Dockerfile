# Use R base image with Bioconductor
FROM rocker/r-ver:4.3.2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('SummarizedExperiment', 'S4Vectors', 'Biobase'), ask=FALSE, update=FALSE)" && \
    R -e "install.packages(c('Matrix', 'jsonlite'), repos='https://cloud.r-project.org')"

# Set up Python virtual environment
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install Python dependencies
WORKDIR /app
COPY server/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY server/ ./server/

# Set working directory to server
WORKDIR /app/server

# App Runner / ECS expects PORT environment variable
ENV PORT=8080

# Run the application
CMD ["python", "main.py"]
