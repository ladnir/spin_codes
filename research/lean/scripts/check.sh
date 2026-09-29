#!/usr/bin/env bash
# Loop invariant check. Must pass at the end of every iteration.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
fail=0

echo "== 1. build (includes SpinCodes.Pin: statement drift is a build error) =="
if ! lake build 2>&1 | tail -1 | grep -q "Build completed successfully"; then
  echo "FAIL: build"; fail=1
else
  echo "ok"
fi

echo "== 2. no sorry / admit / project axiom / native_decide =="
# native_decide would add Lean.ofReduceBool to the axiom closure: the
# numerical certificates must be kernel-checked, not compiled-evaluated.
if grep -rnE '\bsorry\b|\badmit\b|^axiom |^ *axiom |\bnative_decide\b' --include=*.lean SpinCodes/ SpinCodes.lean; then
  echo "FAIL: sorry/admit/axiom/native_decide present"; fail=1
else
  echo "ok"
fi

echo "== 3. axiom closure of the headline theorem =="
tmp="$(mktemp -t axcheck.XXXXXX)".lean
printf 'import SpinCodes\n#print axioms Spin.Structured.ScalableCertificate.distance_whp\n' > "$tmp"
out="$(lake env lean "$tmp" 2>&1)"
rm -f "$tmp"
echo "$out"
if ! echo "$out" | grep -q 'depends on axioms: \[propext, Classical.choice, Quot.sound\]$'; then
  echo "FAIL: unexpected axiom closure"; fail=1
else
  echo "ok"
fi

echo "== 4. numerical certificates are kernel-checked =="
tmp2="$(mktemp -t certcheck.XXXXXX)".lean
printf 'import SpinCodes
#print axioms Spin.Structured.sparseRows_cert
' > "$tmp2"
out2="$(lake env lean "$tmp2" 2>&1)"
rm -f "$tmp2"
echo "$out2"
if echo "$out2" | grep -q 'ofReduceBool\|trustCompiler'; then
  echo "FAIL: certificate relies on compiled evaluation"; fail=1
else
  echo "ok"
fi

if [ "$fail" -eq 0 ]; then echo "ALL CHECKS PASS"; else echo "CHECKS FAILED"; fi
exit "$fail"
