#!/usr/bin/env bash
# Replay the complete T3a Collatz bridge and its repeated-step bound.
set -euo pipefail
cd "$(dirname "$0")/.."
python -X utf8 -B scripts/sparse_bridge_data.py
python -X utf8 -B scripts/sparse_dominance_data.py
python -X utf8 -B scripts/sparse_maxima.py
python -X utf8 -B scripts/sparse_residual_data.py
python -X utf8 -B scripts/sparse_residual_sound.py
# Finish dependency builds before the numerical workers start.
lake build SpinCodes SpinCodes.Structured.PolyPacked SpinCodes.Structured.SparseSelection \
  SpinCodes.Structured.SparseWeights SpinCodes.Structured.SparseProgramDefs
python -X utf8 -B scripts/check-sparse-bridge.py --jobs 4 --phase weights
python -X utf8 -B scripts/check-sparse-bridge.py --jobs 4 --phase dominance
python -X utf8 -B scripts/check-sparse-bridge.py --jobs 4 --phase contributions
for module in SparseMaxima SparseMaximumAll SparseContributionSound SparseResiduals \
  SparseResidualsSound SparseProgram; do
  lake env lean -o ".lake/build/lib/lean/SpinCodes/Structured/$module.olean" \
    "SpinCodes/Structured/$module.lean"
done
python -X utf8 -B scripts/check-sparse-semantic.py
lake env lean scripts/sparse_axioms.lean
