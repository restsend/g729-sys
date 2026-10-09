#!/bin/bash
# Fuzz testing runner for g729-sys.
# Usage: ./fuzz.sh [target] [duration_seconds] [extra cargo-fuzz args...]
set -e

TARGET=${1:-fuzz_decoder}
DURATION=${2:-60}
shift 2 2>/dev/null || true

if ! cargo fuzz --version >/dev/null 2>&1; then
    echo "Installing cargo-fuzz..."
    cargo install cargo-fuzz
fi

cd "$(dirname "$0")/.."
cargo +nightly fuzz run "$TARGET" --release -- -max_total_time="$DURATION" "$@"
