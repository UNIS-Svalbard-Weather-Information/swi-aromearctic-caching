# Stage 1: Dependency resolution
FROM astral/uv:python3.13-bookworm-slim AS uv
WORKDIR /swi
COPY pyproject.toml uv.lock ./
# Export the locked, hashed set: the image installs exactly what uv.lock pins.
RUN uv export --frozen --no-dev --no-emit-project -o requirements.txt

# Stage 2: Build stage
FROM python:3.13-slim AS builder
WORKDIR /swi
COPY --from=uv /swi/requirements.txt .
RUN pip install --no-cache-dir --require-hashes -r requirements.txt

# Stage 3: Runtime stage
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
# The job never installs packages at runtime: remove pip (and its vendored copies of
# urllib3, msgpack, ...) so they are neither exposed nor flagged by image scanners.
RUN python -m pip uninstall -y pip setuptools wheel 2>/dev/null || true
COPY . .

# Make the entrypoint scripts executable
RUN chmod +x /swi/entrypoint.sh /swi/run-cron.sh /swi/main.py

# Use the entrypoint script - supports DOCKER_CRON=1 to stay alive between runs so an
# external scheduler can re-trigger the job via `docker exec <container> run-cron`
ENTRYPOINT ["./entrypoint.sh"]
