#!/bin/bash
# Wrapper script for running BAMCP MCP server with virtualenv
exec /home/rtj/workspace/BAMCP/.venv/bin/python -m bamcp "$@"
