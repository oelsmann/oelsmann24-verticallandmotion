FROM python:3.12-slim-bookworm
COPY --from=ghcr.io/astral-sh/uv:0.7.11 /uv /uvx /bin/

# Where we're installing this thing.
ARG APP_HOME="/opt/oelsmann24-verticallandmotion"

WORKDIR ${APP_HOME}
COPY . .

# Install the application dependencies into a local virtual environment, compiling to bytecode.
RUN uv sync --frozen --no-cache --no-dev --compile-bytecode

# Easily run commands from the environment just created.
ENV PATH="${APP_HOME}/.venv/bin:$PATH"

# Use custom user/group so the application does not run with root permissions.
USER 9876:9876

ENTRYPOINT ["oelsmann24-verticallandmotion"]
