# Stage 1: Build stage
FROM python:3.14-slim AS builder

WORKDIR /swi

# Install uv
RUN pip install uv

# Copy only the necessary files for dependency installation
COPY pyproject.toml .
COPY uv.lock .

# Install dependencies in a virtual environment
RUN uv pip install --system --no-cache-dir .

# Stage 2: Runtime stage
FROM python:3.14-slim

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

# Make the entrypoint scripts executable
RUN chmod +x /swi/entrypoint.sh /swi/run-cron.sh /swi/main.py

# Use the entrypoint script - supports DOCKER_CRON=1 to stay alive between runs so an
# external scheduler can re-trigger the job via `docker exec <container> run-cron`
ENTRYPOINT ["./entrypoint.sh"]
