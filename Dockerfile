# Stage 1: Build stage
FROM python:3.13-slim AS builder

WORKDIR /swi

# Install uv
RUN pip install uv

# Copy only the necessary files for dependency installation
COPY pyproject.toml .
COPY uv.lock .

# Install dependencies in a virtual environment
RUN uv pip install --system --no-cache-dir .

# Stage 2: Runtime stage
FROM python:3.13-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libgdal-dev \
    libproj-dev \
    && rm -rf /var/lib/apt/lists/*


WORKDIR /swi

# Copy only the installed dependencies and your application code
COPY --from=builder /usr/local/lib/python3.13/site-packages /usr/local/lib/python3.13/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin
COPY . .

# Ensure the script is executable
RUN chmod +x /swi/main.py

# Run the script
CMD ["python", "/swi/main.py"]
