#!/usr/bin/env bash
# Build and run the OT-generation kernels on a Zen 4 host, serially, pinned to
# one core, under the shared benchmark locks, and only when no build runs.
# Usage: run_peach.sh <permute_conv checkout> [core] [reps]
set -euo pipefail
ROOT=${1:?checkout}; CORE=${2:-15}; REPS=${3:-21}
B=${OTGEN_BUILD:-$HOME/turbo/otgen-build}
busy() { pgrep -a -x 'lean|lake|cc1plus|cargo|rustc|ninja|make' || true; pgrep -fa 'bench|frontend_libOTe|taskset' | grep -v pgrep || true; }
cmake -S "$ROOT/research/workstreams/otgen_kernels" -B "$B" -DCMAKE_BUILD_TYPE=Release "-DOTGEN_ARCH=-march=znver4" >/dev/null
cmake --build "$B" -j8 2>&1 | grep -v ignored-attributes | grep -E "error|FAILED" || true
test -x "$B/otgen_bench" && test -x "$B/otgen_fused" && test -x "$B/otgen_fused_outline"
for i in $(seq 1 120); do [ -z "$(busy)" ] && break; sleep 30; done
[ -z "$(busy)" ] || { echo "host busy; not benchmarking"; busy | cut -c1-120; exit 1; }
LK="flock /tmp/prindal-addition-encoder-benchmark.lock flock /tmp/bare-spin-benchmark.lock flock /tmp/hypercat-benchmark.lock flock /tmp/hypercat-global-benchmark.lock"
if [ "${OTGEN_SKIP_BENCH:-0}" != 1 ]; then echo "== otgen_bench"; $LK taskset -c "$CORE" "$B/otgen_bench" "$REPS"; fi
echo "== otgen_fused"; $LK taskset -c "$CORE" "$B/otgen_fused" "$REPS"
[ "${OTGEN_SKIP_OUTLINE:-0}" = 1 ] || { echo "== otgen_fused_outline"; $LK taskset -c "$CORE" "$B/otgen_fused_outline" "$REPS"; }
