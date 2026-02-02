#!/bin/sh
# BAMCP Docker entrypoint
# Runs the health check first, then executes the MCP server (or a custom command).

set -e

# If the first argument is a flag (starts with -), prepend the default command.
if [ "${1#-}" != "$1" ]; then
    set -- python -m bamcp "$@"
fi

# If the command is the default server, run a quick pre-flight check.
if [ "$1" = "python" ] && [ "$2" = "-m" ] && [ "$3" = "bamcp" ]; then
    echo "Running pre-flight health check..." >&2
    python /app/docker/healthcheck.py || {
        echo "Health check failed, aborting." >&2
        exit 1
    }
    echo "Health check passed." >&2

    transport="${BAMCP_TRANSPORT:-stdio}"
    if [ "$transport" != "stdio" ]; then
        echo "Starting BAMCP (transport=${transport}, port=${BAMCP_PORT:-8000})..." >&2
    else
        echo "Starting BAMCP (transport=stdio)..." >&2
    fi
fi

exec "$@"
