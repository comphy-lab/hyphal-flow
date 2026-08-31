#!/bin/bash
# Thin wrapper. Prefer runSimulation.sh from the repository root.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NP="${2:-4}"

echo "Use runSimulation.sh instead of this wrapper."
echo "Example: bash ${ROOT}/runSimulation.sh default.params --mpi --CPUs ${NP}"
exit 1
