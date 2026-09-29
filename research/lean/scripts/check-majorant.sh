#!/usr/bin/env bash
# Replay the full refined BA majorant, including its real theorem and pins.
# Run from research/lean. Generated data must already be present.
set -euo pipefail
cd "$(dirname "$0")/.."
lake build SpinCodes
python -X utf8 -B scripts/check-majorant.py --phase parts --jobs "${MAJORANT_JOBS:-6}"
python -X utf8 -B scripts/check-majorant.py --phase assembly --jobs 2
for name in RefinedData Refined Pin; do
  lake env lean -o ".lake/build/lib/lean/SpinCodes/Majorant/$name.olean" "SpinCodes/Majorant/$name.lean"
done
out=$(lake env lean scripts/majorant_axioms.lean)
printf '%s\n' "$out"
if echo "$out" | grep -qE 'sorryAx|ofReduceBool|trustCompiler'; then
  echo 'FAIL: forbidden axiom in majorant proof'
  exit 1
fi
echo "$out" | grep -q 'Spin.Majorant.baExponent_le_refined.*depends on axioms:'
echo 'MAJORANT CHECKS PASS'
