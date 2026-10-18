# Stage 1: Dependency resolution
FROM astral/uv:python3.13-bookworm-slim AS uv
WORKDIR /swi
COPY pyproject.toml uv.lock ./
# Export the locked, hashed set: the image installs exactly what uv.lock pins.
RUN uv export --frozen --no-dev --no-emit-project -o requirements.txt

# Stage 2: Build stage
FROM python:3.14-slim AS builder
WORKDIR /swi
COPY --from=uv /swi/requirements.txt .
RUN pip install --no-cache-dir --require-hashes -r requirements.txt

# Stage 3: Runtime stage
FROM python:3.14-slim
# The project reserves uids above 9000. Every SWI image that shares a volume with this one
# (the other cron jobs, mapproxy, the read-only data services) uses the same uid and gid,
# so files written here are owned by the same user everywhere.
ARG APP_UID=9001
ARG APP_GID=9001

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

# Unprivileged user with a home directory (entrypoint.sh writes ~/.bashrc)
RUN groupadd --gid "${APP_GID}" swi \
    && useradd --uid "${APP_UID}" --gid "${APP_GID}" --create-home --shell /usr/sbin/nologin swi

# The job writes next to its code (app.log, ./temp) and into the forecast volume. Create the
# mount point owned by the job user: an empty named volume inherits the ownership of the
# directory it is mounted on, so a fresh volume is writable without any runtime chown.
RUN mkdir -p /swi/data/forecast \
    && chown -R swi:swi /swi

# Copy application code (stays owned by root: the job can create files here, not edit code)
COPY . .

# Make the entrypoint scripts executable
RUN chmod +x /swi/entrypoint.sh /swi/run-cron.sh /swi/main.py

USER swi

# Use the entrypoint script - supports DOCKER_CRON=1 to stay alive between runs so an
# external scheduler can re-trigger the job via `docker exec <container> run-cron`
ENTRYPOINT ["./entrypoint.sh"]
