#!/usr/bin/env python3
"""Exact-local drive-stratified masses and point caps for packet-8 Riffle.

For packet-class fugacities f_j and an incoming state support R of weight q,
write U for the current 64-bit input and W=U xor R for the accumulator drive.
This probe computes

    H[q,d,y] = sum f(U)

exactly over pairs with wt(R)=q, wt(W)=d, and wt(Acc(W))=y.  It also computes
the exact maximum, over drives W of weight d, of

    c[q,d] = C(64,q)^-1 sum_{R:wt(R)=q} f(W xor R).

The maximization is exact because each byte polynomial depends only on the
drive-byte weight; all 12870 multisets of eight byte weights are enumerated.
Binary64 fugacities make this a diagnostic table.  The same finite formulas
can later be evaluated with rational/outward arithmetic for a certificate.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from probe_packet8_hard_profile_thermodynamic import PACKET_WIDTH
from probe_packet8_weight_profile_scalar import CLASSES


BITS = 64
PACKETS = BITS // PACKET_WIDTH


def one_packet_transfer(fugacities: np.ndarray) -> np.ndarray:
    transfer = np.zeros((2, 2, 9, 9, 9), dtype=np.float64)
    for current_input in range(256):
        coefficient = float(fugacities[current_input.bit_count()])
        if coefficient == 0.0:
            continue
        for selected_state in range(256):
            selected = selected_state.bit_count()
            drive = current_input ^ selected_state
            drive_weight = drive.bit_count()
            for incoming_parity in range(2):
                parity = incoming_parity
                emitted = 0
                for position in range(PACKET_WIDTH):
                    parity ^= (drive >> position) & 1
                    emitted += parity
                transfer[
                    incoming_parity,
                    parity,
                    selected,
                    drive_weight,
                    emitted,
                ] += coefficient
    return transfer


def _up_add(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Componentwise IEEE-754 upper addition for nonnegative arrays.

    Exact zeros are preserved so outward evaluation does not create spurious
    subnormal support.  Every positive rounded result is advanced by one ulp,
    which encloses the exact sum of the represented binary rationals.
    """

    raw = left + right
    return np.where(right > 0.0, np.nextafter(raw, np.inf), left)


def _up_mul(left: np.ndarray | float, right: np.ndarray | float) -> np.ndarray:
    left_array = np.asarray(left)
    right_array = np.asarray(right)
    raw = left_array * right_array
    positive_underflow = (left_array > 0.0) & (right_array > 0.0) & (raw == 0.0)
    if np.any(positive_underflow):
        raise FloatingPointError(
            "outward multiplication underflowed a structurally positive product"
        )
    return np.where(raw > 0.0, np.nextafter(raw, np.inf), raw)


def one_packet_transfer_outward(fugacities: np.ndarray) -> np.ndarray:
    """Outward upper version of :func:`one_packet_transfer`.

    The supplied binary64 fugacities are treated as exact dyadic rationals.
    Only nonnegative additions occur in this first layer.
    """

    transfer = np.zeros((2, 2, 9, 9, 9), dtype=np.float64)
    for current_input in range(256):
        coefficient = float(fugacities[current_input.bit_count()])
        if coefficient == 0.0:
            continue
        for selected_state in range(256):
            selected = selected_state.bit_count()
            drive = current_input ^ selected_state
            drive_weight = drive.bit_count()
            for incoming_parity in range(2):
                parity = incoming_parity
                emitted = 0
                for position in range(PACKET_WIDTH):
                    parity ^= (drive >> position) & 1
                    emitted += parity
                index = (
                    incoming_parity,
                    parity,
                    selected,
                    drive_weight,
                    emitted,
                )
                transfer[index] = np.nextafter(
                    transfer[index] + coefficient, np.inf
                )
    return transfer


def block_histograms(fugacities: np.ndarray) -> np.ndarray:
    packet = one_packet_transfer(fugacities)
    entries = [
        (incoming, outgoing, selected, drive, emitted, packet[incoming, outgoing, selected, drive, emitted])
        for incoming in range(2)
        for outgoing in range(2)
        for selected in range(9)
        for drive in range(9)
        for emitted in range(9)
        if packet[incoming, outgoing, selected, drive, emitted]
    ]
    dp = np.zeros((2, BITS + 1, BITS + 1, BITS + 1), dtype=np.float64)
    dp[0, 0, 0, 0] = 1.0
    maximum = 0
    for _ in range(PACKETS):
        next_dp = np.zeros_like(dp)
        for incoming, outgoing, selected, drive, emitted, coefficient in entries:
            source = dp[incoming, : maximum + 1, : maximum + 1, : maximum + 1]
            next_dp[
                outgoing,
                selected : selected + maximum + 1,
                drive : drive + maximum + 1,
                emitted : emitted + maximum + 1,
            ] += coefficient * source
        dp = next_dp
        maximum += PACKET_WIDTH
    result = dp.sum(axis=0)
    packet_mass = sum(
        classes * float(fugacity)
        for classes, fugacity in zip(CLASSES, fugacities)
    )
    for state_weight in range(BITS + 1):
        expected = math.comb(BITS, state_weight) * packet_mass**PACKETS
        actual = float(np.sum(result[state_weight]))
        if abs(actual - expected) > 3e-12 * max(1.0, expected):
            raise SystemExit(
                f"drive-stratified cap: mass mismatch at q={state_weight}"
            )
    return result


def block_histograms_outward(fugacities: np.ndarray) -> np.ndarray:
    """Componentwise upper enclosure of the exact local histogram."""

    packet = one_packet_transfer_outward(fugacities)
    entries = [
        (incoming, outgoing, selected, drive, emitted, packet[incoming, outgoing, selected, drive, emitted])
        for incoming in range(2)
        for outgoing in range(2)
        for selected in range(9)
        for drive in range(9)
        for emitted in range(9)
        if packet[incoming, outgoing, selected, drive, emitted]
    ]
    dp = np.zeros((2, BITS + 1, BITS + 1, BITS + 1), dtype=np.float64)
    dp[0, 0, 0, 0] = 1.0
    maximum = 0
    for _ in range(PACKETS):
        next_dp = np.zeros_like(dp)
        for incoming, outgoing, selected, drive, emitted, coefficient in entries:
            source = dp[incoming, : maximum + 1, : maximum + 1, : maximum + 1]
            contribution = _up_mul(coefficient, source)
            destination = next_dp[
                outgoing,
                selected : selected + maximum + 1,
                drive : drive + maximum + 1,
                emitted : emitted + maximum + 1,
            ]
            destination[...] = _up_add(destination, contribution)
        dp = next_dp
        maximum += PACKET_WIDTH
    # The two terminal accumulator parities are disjoint alternatives.
    return _up_add(dp[0], dp[1])


def byte_polynomials(fugacities: np.ndarray) -> tuple[np.ndarray, ...]:
    rows = []
    for drive_weight in range(9):
        polynomial = np.zeros(9, dtype=np.float64)
        for selected in range(9):
            for overlap in range(
                max(0, selected - (PACKET_WIDTH - drive_weight)),
                min(drive_weight, selected) + 1,
            ):
                current_weight = drive_weight + selected - 2 * overlap
                count = math.comb(drive_weight, overlap) * math.comb(
                    PACKET_WIDTH - drive_weight, selected - overlap
                )
                polynomial[selected] += count * float(fugacities[current_weight])
        rows.append(polynomial)
    return tuple(rows)


def exact_point_caps(fugacities: np.ndarray) -> np.ndarray:
    local = byte_polynomials(fugacities)
    maxima = np.zeros((BITS + 1, BITS + 1), dtype=np.float64)
    for weights in itertools.combinations_with_replacement(range(9), PACKETS):
        polynomial = np.array([1.0])
        for weight in weights:
            polynomial = np.convolve(polynomial, local[weight])
        drive_weight = sum(weights)
        maxima[:, drive_weight] = np.maximum(maxima[:, drive_weight], polynomial)
    for state_weight in range(BITS + 1):
        maxima[state_weight] /= math.comb(BITS, state_weight)
    return maxima


def byte_polynomials_outward(fugacities: np.ndarray) -> tuple[np.ndarray, ...]:
    rows = []
    for drive_weight in range(9):
        polynomial = np.zeros(9, dtype=np.float64)
        for selected in range(9):
            for overlap in range(
                max(0, selected - (PACKET_WIDTH - drive_weight)),
                min(drive_weight, selected) + 1,
            ):
                current_weight = drive_weight + selected - 2 * overlap
                count = math.comb(drive_weight, overlap) * math.comb(
                    PACKET_WIDTH - drive_weight, selected - overlap
                )
                term = _up_mul(float(count), fugacities[current_weight])
                polynomial[selected] = np.nextafter(
                    polynomial[selected] + float(term), np.inf
                )
        rows.append(polynomial)
    return tuple(rows)


def _convolve_outward(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.zeros(len(left) + len(right) - 1, dtype=np.float64)
    for left_index, left_value in enumerate(left):
        if left_value == 0.0:
            continue
        contribution = _up_mul(left_value, right)
        target = result[left_index : left_index + len(right)]
        target[...] = _up_add(target, contribution)
    return result


def exact_point_caps_outward(fugacities: np.ndarray) -> np.ndarray:
    """Componentwise upper enclosure of every exact drive point cap."""

    local = byte_polynomials_outward(fugacities)
    maxima = np.zeros((BITS + 1, BITS + 1), dtype=np.float64)
    for weights in itertools.combinations_with_replacement(range(9), PACKETS):
        polynomial = np.array([1.0])
        for weight in weights:
            polynomial = _convolve_outward(polynomial, local[weight])
        drive_weight = sum(weights)
        maxima[:, drive_weight] = np.maximum(maxima[:, drive_weight], polynomial)
    for state_weight in range(BITS + 1):
        positive = maxima[state_weight] > 0.0
        quotient = maxima[state_weight, positive] / math.comb(BITS, state_weight)
        maxima[state_weight, positive] = np.nextafter(quotient, np.inf)
    return maxima


def empirical_cap_summary(trace_path: Path, caps: np.ndarray, global_cap: float) -> dict:
    data = json.loads(trace_path.read_text(encoding="utf-8"))
    traces = data.get("transition_traces", [])
    if not traces:
        raise SystemExit("drive-stratified cap: trace has no transition data")
    trace = max(traces, key=lambda row: float(row["beta"]))
    savings: list[tuple[float, int]] = []
    missing = 0
    for row in trace["drive_transition"]:
        q = int(row["state_weight"])
        drive = int(row["drive_weight"])
        count = int(row["count"])
        cap = float(caps[q, drive])
        if cap <= 0.0:
            missing += count
            continue
        savings.append((math.log2(global_cap / cap), count))
    total = sum(count for _value, count in savings)
    ordered = sorted(savings)

    def quantile(fraction: float) -> float:
        target = fraction * total
        cumulative = 0
        for value, count in ordered:
            cumulative += count
            if cumulative >= target:
                return value
        return ordered[-1][0]

    return {
        "source": str(trace_path),
        "beta": trace["beta"],
        "observations": int(trace["observations"]),
        "covered_observations": total,
        "missing_observations": missing,
        "mean_log2_cap_saving": sum(value * count for value, count in savings) / total,
        "saving_quantiles": {
            "q10": quantile(0.10),
            "q25": quantile(0.25),
            "q50": quantile(0.50),
            "q75": quantile(0.75),
            "q90": quantile(0.90),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fugacities", required=True)
    parser.add_argument("--trace", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    fugacities = np.array(
        [float(value) for value in args.fugacities.split(",")], dtype=np.float64
    )
    if len(fugacities) != 9 or np.any(fugacities < 0.0) or fugacities[0] <= 0.0:
        raise SystemExit("drive-stratified cap: invalid fugacities")

    histograms = block_histograms(fugacities)
    caps = exact_point_caps(fugacities)
    global_cap = float(np.max(fugacities)) ** PACKETS
    positive = caps[caps > 0.0]
    payload: dict[str, object] = {
        "fugacities": fugacities.tolist(),
        "global_point_cap": global_cap,
        "positive_cap_min": float(np.min(positive)),
        "positive_cap_max": float(np.max(positive)),
        "histogram_shape": list(histograms.shape),
        "cap_shape": list(caps.shape),
    }
    if args.trace:
        payload["empirical"] = empirical_cap_summary(args.trace, caps, global_cap)

    print("packet-8 exact-local drive-stratified cap probe")
    print("histogram_axes=incoming_state_weight,drive_weight,emitted_weight")
    print(f"histogram_shape={histograms.shape}")
    print(f"global_point_cap_log2={math.log2(global_cap):.12f}")
    print(
        f"positive_stratified_cap_log2_range="
        f"{math.log2(float(np.min(positive))):.12f},"
        f"{math.log2(float(np.max(positive))):.12f}"
    )
    if args.trace:
        empirical = payload["empirical"]
        print(
            f"trace_beta={empirical['beta']:.9f} "
            f"mean_log2_cap_saving={empirical['mean_log2_cap_saving']:.9f} "
            f"saving_quantiles="
            + ",".join(
                f"{name}:{value:.6f}"
                for name, value in empirical["saving_quantiles"].items()
            )
        )
    print("status=EXACT_LOCAL_COMBINATORICS_BINARY64_FUGACITY_DIAGNOSTIC")

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.output.with_suffix(".npz"), histograms=histograms, caps=caps
        )
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"output={args.output}")
        print(f"table_output={args.output.with_suffix('.npz')}")


if __name__ == "__main__":
    main()
