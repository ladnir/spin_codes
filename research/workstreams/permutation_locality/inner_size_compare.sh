#!/usr/bin/env bash
# Usage: bash inner_size_compare.sh ROOT EXPONENT VARIANT:MODE VARIANT:MODE [...]
# base selects ROOT/inner-size; NAME selects ROOT/inner-size-NAME.
# Optional environment: SPIN_COMPARE_CALLS (301), SPIN_COMPARE_SEEDS
# (space-separated: 1 17 43 91), SPIN_COMPARE_PROFILE (0 or 1, default 0).
set -euo pipefail

root=$(realpath -- "${1:?comparison directory}")
exponent=${2:?message exponent}
shift 2
case "$exponent" in
  14|15|16|17|18|19|20) ;;
  *) printf 'Exponent must be in 14..20\n' >&2; exit 2 ;;
esac
if (( $# < 2 )); then
  printf 'At least two explicit VARIANT:MODE selections are required\n' >&2
  exit 2
fi

calls=${SPIN_COMPARE_CALLS-301}
profile=${SPIN_COMPARE_PROFILE-0}
seed_text=${SPIN_COMPARE_SEEDS-'1 17 43 91'}
if [[ ! $calls =~ ^[1-9][0-9]*$ ]]; then
  printf 'SPIN_COMPARE_CALLS must be a positive decimal integer\n' >&2
  exit 2
fi
if [[ $profile != 0 && $profile != 1 ]]; then
  printf 'SPIN_COMPARE_PROFILE must be 0 or 1\n' >&2
  exit 2
fi
seed_text=${seed_text//$'\n'/ }
read -r -a seeds <<< "$seed_text"
if (( ${#seeds[@]} == 0 )); then
  printf 'SPIN_COMPARE_SEEDS must contain at least one seed\n' >&2
  exit 2
fi
declare -A seen_seeds=()
for seed in "${seeds[@]}"; do
  if [[ ! $seed =~ ^(0|[1-9][0-9]*)$ || -n ${seen_seeds[$seed]+present} ]]; then
    printf 'Seeds must be distinct canonical nonnegative decimal integers: %s\n' "$seed" >&2
    exit 2
  fi
  seen_seeds[$seed]=1
done

variants=()
modes=()
binaries=()
declare -A seen_selections=()
for selection in "$@"; do
  if [[ ! $selection =~ ^([A-Za-z0-9-]+):(0|[1-9][0-9]*)$ ]]; then
    printf 'Invalid VARIANT:MODE selection: %s\n' "$selection" >&2
    exit 2
  fi
  variant=${BASH_REMATCH[1]}
  mode=${BASH_REMATCH[2]}
  if [[ -n ${seen_selections[$selection]+present} ]]; then
    printf 'Duplicate selection would overwrite its logs: %s\n' "$selection" >&2
    exit 2
  fi
  seen_selections[$selection]=1
  binary="$root/inner-size-$variant"
  if [[ $variant == base ]]; then
    binary="$root/inner-size"
  fi
  if [[ ! -x $binary ]]; then
    printf 'Missing executable: %s\n' "$binary" >&2
    exit 2
  fi
  variants+=("$variant")
  modes+=("$mode")
  binaries+=("$binary")
done

exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -n 7 || exit 75

logs=$(mktemp -d "$root/compare-XXXXXX")
printf 'Logs: %s\n' "$logs"
printf 'exponent=%s\ncalls=%s\nprofile=%s\nseeds=%s\nrepeats=2\ncpu=15\n' \
  "$exponent" "$calls" "$profile" "${seeds[*]}" > "$logs/context.txt"
printf 'forward_order=%s\n' "$*" >> "$logs/context.txt"
printf 'repeat_0=forward\nrepeat_1=reverse\n' >> "$logs/context.txt"
sha256sum "$0" > "$logs/hashes.txt"
declare -A hashed_binaries=()
for binary in "${binaries[@]}"; do
  if [[ -z ${hashed_binaries[$binary]+present} ]]; then
    sha256sum "$binary" >> "$logs/hashes.txt"
    hashed_binaries[$binary]=1
  fi
done

for seed in "${seeds[@]}"; do
  for ((rep=0;rep<2;++rep)); do
    for ((slot=0;slot<${#variants[@]};++slot)); do
      index=$slot
      if ((rep)); then
        index=$((${#variants[@]}-1-slot))
      fi
      variant=${variants[$index]}
      mode=${modes[$index]}
      binary=${binaries[$index]}
      output="$logs/full-$variant-$mode-$exponent-$seed-$rep.csv"
      taskset -c 15 "$binary" "$mode" "$exponent" "$seed" "$calls" "$profile" > "$output"
      cat "$output"
    done
  done
done
