#!/bin/sh
# BAMCP Docker entrypoint
# Runs the health check first, then executes the MCP server (or a custom command).

set -e

# If the first argument is a flag (starts with -), prepend the default command.
if [ "${1#-}" != "$1" ]; then
    set -- python -m bamcp "$@"
fi

# If the command is the default server, run a quick health check first.
if [ "$1" = "python" ] && [ "$2" = "-m" ] && [ "$3" = "bamcp" ]; then
    echo "Running pre-flight health check..." >&2
    python /app/docker/healthcheck.py || {
        echo "Health check failed, aborting." >&2
        exit 1
    }
    echo "Health check passed." >&2
fi

exec "$@"
