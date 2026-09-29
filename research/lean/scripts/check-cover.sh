#!/usr/bin/env bash
# Check the interval covers.
#
# These are deliberately NOT part of `scripts/check.sh`.  The small cover takes
# about fourteen minutes and the full one about an hour and a half, which is
# not a cost worth paying on every build — but it does mean a cover's claim is
# certified only when this script has passed.  Record each passing run in
# LOOP_LOG.md with its date and the module checked.
#
# Usage:  bash scripts/check-cover.sh [small|full|all]

set -euo pipefail
cd "$(dirname "$0")/.."

WHICH="${1:-all}"

check_one () {
  local mod="$1"
  echo "== building $mod =="
  local t0 t1
  t0=$(date +%s)
  lake build "$mod"
  t1=$(date +%s)
  echo "ok ($((t1 - t0))s)"
}

axioms_of () {
  local mod="$1" name="$2"
  local tmp
  tmp=$(mktemp --suffix=.lean)
  printf 'import %s\n#print axioms %s\n' "$mod" "$name" > "$tmp"
  local out
  out=$(lake env lean "$tmp")
  rm -f "$tmp"
  echo "$out"
  if echo "$out" | grep -qE "ofReduceBool|trustCompiler|sorryAx"; then
    echo "FAILED: $name has a forbidden axiom"
    exit 1
  fi
  if echo "$out" | grep -qE "propext|Classical.choice|Quot.sound|does not depend"; then
    echo "ok"
  else
    echo "FAILED: unexpected axiom line for $name"
    exit 1
  fi
}

if [ "$WHICH" = "small" ] || [ "$WHICH" = "all" ]; then
  check_one SpinCodes.Cover.DenseTailSmall
  axioms_of SpinCodes.Cover.DenseTailSmall Spin.Cover.denseTailSmall
fi

if [ "$WHICH" = "full" ] || [ "$WHICH" = "all" ]; then
  if [ -f SpinCodes/Cover/DenseTailReal.lean ]; then
    # The 44 parts are independent; build them in batches so at most a few
    # Lean processes run at once.  Lake has no job-count flag in this version
    # and would otherwise use every core, which at a few GB a part exceeds
    # the machine's memory.
    N=$(ls SpinCodes/Cover/Part[0-9]*.lean | wc -l)
    i=0
    while [ "$i" -lt "$N" ]; do
      targets=""
      for k in 0 1 2 3 4 5; do
        j=$((i+k))
        [ "$j" -lt "$N" ] && targets="$targets SpinCodes.Cover.Part$j"
      done
      echo "== parts $i.. =="
      lake build $targets
      i=$((i+6))
    done
    check_one SpinCodes.Cover.DenseTailReal
    axioms_of SpinCodes.Cover.DenseTailReal Spin.Cover.denseTail_real
  else
    echo "SpinCodes/Cover/DenseTailReal.lean not generated; skipping"
  fi
fi

echo "COVER CHECKS PASS"
