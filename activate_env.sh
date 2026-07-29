#!/usr/bin/env bash
# Activate the tagseq Python environment
ENV_DIR="$(cd "$(dirname "$0")" && pwd)"
if [ -d "$ENV_DIR/.venv" ]; then
    source "$ENV_DIR/.venv/bin/activate"
elif [ -d "$HOME/miniconda3/envs/tagseq" ]; then
    source "$HOME/miniconda3/bin/activate" tagseq
else
    echo "No tagseq environment found. Run: uv venv && uv pip install -e ."
fi
export PATH="$ENV_DIR/bin:$PATH"
