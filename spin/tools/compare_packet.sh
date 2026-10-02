#!/usr/bin/env bash
# Serial relative-performance guard for the precomputed packet transpose.
#
# Usage: bash spin/tools/compare_packet.sh REFERENCE CANDIDATE EXISTING_LOG_DIR
# Both executables must implement the existing spin_packet_bench CLI/output.
# Build the same bench/packet.cpp against the two libraries, using matched
# compiler, Release, ISA, and SPIN_TUNE settings. Research harnesses with a
# different CLI or warmup count are not interchangeable with this protocol.
# No executable is built here; setup/allocation remain outside its timed calls.
#
# Every run covers K=2^16, 2^18, 2^19, and 2^20. For each seed it runs both
# reference/candidate and candidate/reference orders on one pinned CPU, while
# holding all three shared benchmark locks. There are no concurrent children.
# Results are median-of-process-medians, not the fastest observed samples.
# A candidate/reference ratio above the configured limit fails the comparison.
# This is a manual matched-machine guard, not a fixed-time CI requirement.
#
# Environment:
#   SPIN_COMPARE_CPU                 15
#   SPIN_COMPARE_CALLS               301 (1..1000000, same for both binaries)
#   SPIN_COMPARE_SEEDS               '1 17 43 91' (distinct uint64 seeds)
#   SPIN_COMPARE_MEMORY              normal (or huge, same for both binaries)
#   SPIN_COMPARE_MAX_RATIO           1.05 (positive decimal, default 5% slack)
#   SPIN_COMPARE_REFERENCE_BUILD     optional compiler/build-setting notes
#   SPIN_COMPARE_CANDIDATE_BUILD     optional compiler/build-setting notes
# Run normal and huge as separate complete campaigns if both are relevant.
# Raw logs are written only beneath a fresh mktemp directory at the supplied
# destination; choose an ignored or external path, never commit raw results.
set -euo pipefail
export LC_ALL=C

usage() {
    printf 'Usage: bash %s REFERENCE_SPIN_PACKET_BENCH CANDIDATE_SPIN_PACKET_BENCH EXISTING_LOG_DIR\n' "$0"
}
if [[ $# == 1 && $1 == --help ]]; then
    usage
    sed -n '2,/^set -euo pipefail/{ /^#/p; }' "$0"
    exit 0
fi
if (( $# != 3 )); then usage >&2; exit 2; fi
for tool in realpath flock taskset mktemp sha256sum awk tee; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        printf 'Required Linux tool unavailable: %s\n' "$tool" >&2
        exit 2
    fi
done
reference=$(realpath -- "$1")
candidate=$(realpath -- "$2")
log_parent=$(realpath -- "$3")
script=$(realpath -- "$0")
if [[ ! -f $reference || ! -x $reference || ! -f $candidate || ! -x $candidate || ! -d $log_parent ]]; then
    printf 'Two existing executable files and an existing log directory are required\n' >&2
    exit 2
fi
if [[ $reference == "$candidate" ]]; then
    printf 'Reference and candidate must be separate executable paths\n' >&2
    exit 2
fi

cpu=${SPIN_COMPARE_CPU-15}
calls=${SPIN_COMPARE_CALLS-301}
memory=${SPIN_COMPARE_MEMORY-normal}
limit=${SPIN_COMPARE_MAX_RATIO-1.05}
seed_text=${SPIN_COMPARE_SEEDS-'1 17 43 91'}
if [[ ! $cpu =~ ^(0|[1-9][0-9]*)$ || ! $calls =~ ^[1-9][0-9]*$ || ${#calls} -gt 7 ]]; then
    printf 'CPU must be a nonnegative decimal integer; calls must be in 1..1000000\n' >&2
    exit 2
fi
if (( calls > 1000000 )); then
    printf 'SPIN_COMPARE_CALLS must be in 1..1000000\n' >&2
    exit 2
fi
if [[ $memory != normal && $memory != huge ]]; then
    printf 'SPIN_COMPARE_MEMORY must be normal or huge\n' >&2
    exit 2
fi
if [[ ! $limit =~ ^(0|[1-9][0-9]*)(\.[0-9]+)?$ ]] ||
   ! awk -v ratio="$limit" 'BEGIN { exit !(ratio>0 && ratio<=10) }'; then
    printf 'SPIN_COMPARE_MAX_RATIO must be a positive decimal at most 10\n' >&2
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
    if [[ ! $seed =~ ^(0|[1-9][0-9]*)$ || ${#seed} -gt 20 ||
          ( ${#seed} == 20 && $seed > 18446744073709551615 ) ||
          -n ${seen_seeds[$seed]+present} ]]; then
        printf 'Seeds must be distinct canonical uint64 decimal integers: %s\n' "$seed" >&2
        exit 2
    fi
    seen_seeds[$seed]=1
done

exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 9 || { printf 'Another encoder benchmark holds the shared lock\n' >&2; exit 75; }
exec 8>/tmp/bare-spin-benchmark.lock
flock -n 8 || { printf 'Another SPIN benchmark holds the shared lock\n' >&2; exit 75; }
exec 7>/tmp/hypercat-benchmark.lock
flock -n 7 || { printf 'Another Hypercat benchmark holds the shared lock\n' >&2; exit 75; }

logs=$(mktemp -d "$log_parent/packet-compare-XXXXXX")
printf 'Packet comparison logs: %s\n' "$logs"
{
    printf 'reference=%s\ncandidate=%s\n' "$reference" "$candidate"
    printf 'exponents=16 18 19 20\ncpu=%s\ncalls=%s\nmemory=%s\n' "$cpu" "$calls" "$memory"
    printf 'seeds=%s\nrepeat0=reference candidate\nrepeat1=candidate reference\n' "${seeds[*]}"
    printf 'max_candidate_reference_ratio=%s\naggregate=median_of_process_medians\n' "$limit"
    printf 'reference_build=%s\ncandidate_build=%s\n' \
        "${SPIN_COMPARE_REFERENCE_BUILD-unrecorded; supply compiler/Release/ISA/tuning notes}" \
        "${SPIN_COMPARE_CANDIDATE_BUILD-unrecorded; supply compiler/Release/ISA/tuning notes}"
    printf 'expected_benchmark=spin_packet_bench; five warmups; inplace reuse; setup excluded\n'
    uname -a
} > "$logs/context.txt"
sha256sum "$script" "$reference" "$candidate" > "$logs/hashes.txt"
runs="$logs/runs.csv"
printf 'exponent,seed,order,variant,K,calls,backend,memory,median_ms,p10_ms,p90_ms,checksum\n' > "$runs"

for exponent in 16 18 19 20; do
    k=$((1 << exponent))
    for seed in "${seeds[@]}"; do
        for order in 0 1; do
            variants=(reference candidate)
            if (( order )); then variants=(candidate reference); fi
            pair_key=''
            for variant in "${variants[@]}"; do
                binary=$reference
                if [[ $variant == candidate ]]; then binary=$candidate; fi
                result="$logs/$exponent-$seed-$order-$variant.csv"
                printf 'K=2^%s seed=%s order=%s %s\n' "$exponent" "$seed" "$order" "$variant"
                taskset -c "$cpu" "$binary" "$k" "$seed" "$calls" "$memory" \
                    > "$result" 2> "$result.stderr"
                # Strictly require the same existing benchmark contract and
                # requested workload; do not accidentally compare a fallback,
                # different record width, different size, or another harness.
                awk -F, -v exponent="$exponent" -v k="$k" -v seed="$seed" \
                    -v order="$order" -v variant="$variant" -v calls="$calls" -v memory="$memory" '
                    NR==1 {
                        if ($0!="K,seed,calls,backend,memory,median_ms,p10_ms,p90_ms,checksum") bad=1
                        next
                    }
                    NR==2 {
                        if (NF!=9 || "x"$1!="x"k || "x"$2!="x"seed || "x"$3!="x"calls ||
                            $4!="2" || $5!=memory || $6!~/^[0-9]+([.][0-9]+)?$/ ||
                            $7!~/^[0-9]+([.][0-9]+)?$/ || $8!~/^[0-9]+([.][0-9]+)?$/ ||
                            $9!~/^[0-9a-fA-F]+$/ || $6<=0 || $7>$6 || $6>$8) bad=1
                        line=exponent "," seed "," order "," variant "," $1 "," $3 "," $4 "," $5 "," $6 "," $7 "," $8 "," $9
                        next
                    }
                    {bad=1}
                    END {
                        if (NR!=2 || bad) {
                            print "Unexpected packet benchmark CSV/workload/backend: " FILENAME > "/dev/stderr"
                            exit 2
                        }
                        print line
                    }' "$result" >> "$runs"
                key=$(awk -F, 'NR==2 {print $4 "," $9}' "$result")
                if [[ -n $pair_key && $key != "$pair_key" ]]; then
                    printf 'Backend/checksum mismatch for K=2^%s seed=%s order=%s\n' "$exponent" "$seed" "$order" >&2
                    exit 2
                fi
                pair_key=$key
                cat "$result"
            done
        done
    done
done

# Refuse a result collected while somebody rebuilt either executable or edited
# this protocol. All locks remain held through validation and summarization.
sha256sum --check "$logs/hashes.txt" > "$logs/hash-validation.txt"
if awk -F, -v limit="$limit" -v expected="$((${#seeds[@]} * 2))" '
    function median(key, count, sorted, i, j, value) {
        count=counts[key]
        for(i=1;i<=count;i++) {
            value=values[key SUBSEP i]
            for(j=i-1;j>=1 && sorted[j]>value;j--) sorted[j+1]=sorted[j]
            sorted[j+1]=value
        }
        return count%2 ? sorted[(count+1)/2] : (sorted[count/2]+sorted[count/2+1])/2
    }
    NR>1 {key=$1 SUBSEP $4; counts[key]++; values[key SUBSEP counts[key]]=$9+0}
    END {
        print "exponent,K,process_medians_per_binary,reference_ms,candidate_ms,candidate_reference_ratio,max_ratio,status"
        split("16 18 19 20",sizes," ")
        for(i=1;i<=4;i++) {
            e=sizes[i]; a=e SUBSEP "reference"; b=e SUBSEP "candidate"
            if(counts[a]!=expected || counts[b]!=expected) {
                print "Incomplete matched cases at exponent " e > "/dev/stderr"
                exit 2
            }
            reference=median(a); candidate=median(b); ratio=candidate/reference
            status=ratio>limit ? "REGRESSION" : "PASS"
            if(ratio>limit) failed=1
            printf "%d,%.0f,%d,%.6f,%.6f,%.6f,%.6f,%s\n",e,2^e,expected,reference,candidate,ratio,limit,status
        }
        exit failed ? 1 : 0
    }' "$runs" | tee "$logs/summary.csv"; then
    printf 'PASS: all four sizes are within the matched relative-performance limit.\n'
else
    printf 'Comparison failed; inspect %s and repeat noisy campaigns before drawing conclusions.\n' "$logs/summary.csv" >&2
    exit 1
fi
