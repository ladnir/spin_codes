"""Export the coupled q_w,h_w BCH-sandwich model with exact integer coefficients.

This is not a solver.  The JSON is intended as an exact interchange/checking format
for an exact-rational LP/MILP implementation.

Half-spectrum variables:
    q_w = A_w(Q),  h_w = common nonzero Q-coset spectrum,
for even 0 <= w <= 128.

Derived:
    p_w = q_w + 255 h_w
    c_w = q_w + 31 h_w

All constraints are written using integers; senses are eq/ge/le.
"""

from __future__ import annotations

import json
from pathlib import Path
from math import comb

from anchors import L71_HALF, U187_DUAL_HALF
from spectrum_exact import (
    full_from_symmetric_half,
    macwilliams,
    symmetric_kraw_coeff,
    symmetric_moment_coeff,
)

N = 256
KQ = 123
KP = 131
WEIGHTS = list(range(0, 129, 2))
DUAL_WEIGHTS = list(range(0, 129, 2))


def vq(w: int) -> str:
    return f"q_{w}"


def vh(w: int) -> str:
    return f"h_{w}"


def term_dict(**kwargs):
    return {k: int(v) for k, v in kwargs.items() if v}


def add_coeff(d: dict[str, int], name: str, value: int) -> None:
    if value:
        d[name] = d.get(name, 0) + int(value)
        if d[name] == 0:
            del d[name]


def main() -> None:
    L = full_from_symmetric_half(L71_HALF)
    Ld = macwilliams(L, 71)

    Ud = full_from_symmetric_half(U187_DUAL_HALF)
    U = macwilliams(Ud, 69)

    constraints = []
    metadata = {
        "n": N,
        "Q_dimension": KQ,
        "P_dimension": KP,
        "variables": [vq(w) for w in WEIGHTS] + [vh(w) for w in WEIGHTS],
        "derived": {
            "p_w": "q_w + 255 h_w",
            "c_w": "q_w + 31 h_w",
        },
        "notes": [
            "All variables represent half-spectrum coefficients; symmetry supplies 256-w.",
            "Odd primal and dual weights are zero.",
            "Use an exact rational solver for rigorous bounds.",
        ],
    }

    def add(name: str, coeffs: dict[str, int], sense: str, rhs: int):
        constraints.append(
            {"name": name, "coeffs": {k: str(v) for k, v in coeffs.items()},
             "sense": sense, "rhs": str(int(rhs))}
        )

    # Hard support / normalization anchors.
    add("q_0", {vq(0): 1}, "eq", 1)
    add("h_0", {vh(0): 1}, "eq", 0)
    for w in WEIGHTS:
        if 0 < w < 40:
            add(f"q_support_{w}", {vq(w): 1}, "eq", 0)
        if w < 38:
            add(f"h_support_{w}", {vh(w): 1}, "eq", 0)

    # Rigorous affine-orbit lower bounds from Wambach's explicit words.
    add("Wambach_h_38_lower", {vh(38): 1}, "ge", 256)
    add("Wambach_q_40_lower", {vq(40): 1}, "ge", 65280)

    # Stronger rigorous bounds from the calibrated radius-five orbit search.
    # Each retained canonical representative is checked for membership; its
    # full affine orbit is disjoint from every other retained orbit.
    add("orbit_search_h_38_lower", {vh(38): 1}, "ge", 29440)
    add("orbit_search_q_40_lower", {vq(40): 1}, "ge", 12990720)
    add("orbit_search_h_40_lower", {vh(40): 1}, "ge", 163072)
    add("orbit_search_q_42_lower", {vq(42): 1}, "ge", 9204480)
    add("orbit_search_h_42_lower", {vh(42): 1}, "ge", 1067264)
    add("orbit_search_q_44_lower", {vq(44): 1}, "ge", 81991680)
    add("orbit_search_h_44_lower", {vh(44): 1}, "ge", 5977088)
    add("orbit_search_q_46_lower", {vq(46): 1}, "ge", 381757440)
    add("orbit_search_h_46_lower", {vh(46): 1}, "ge", 27937536)
    add("orbit_search_q_48_lower", {vq(48): 1}, "ge", 1980725760)
    add("orbit_search_h_48_lower", {vh(48): 1}, "ge", 114343168)
    add("orbit_search_q_50_lower", {vq(50): 1}, "ge", 679042560)
    add("orbit_search_h_50_lower", {vh(50): 1}, "ge", 28884992)

    # Nonnegativity represented explicitly for portability.
    for w in WEIGHTS:
        add(f"q_nonneg_{w}", {vq(w): 1}, "ge", 0)
        add(f"h_nonneg_{w}", {vh(w): 1}, "ge", 0)

    # OA strength 15 for Q and for every Q coset H.
    for t in range(16):
        cq = {}
        ch = {}
        for w in WEIGHTS:
            a = symmetric_moment_coeff(N, t, w)
            add_coeff(cq, vq(w), a)
            add_coeff(ch, vh(w), a)
        rhs = (1 << (KQ - t)) * comb(N, t)
        add(f"Q_OA_t{t}", cq, "eq", rhs)
        add(f"H_OA_t{t}", ch, "eq", rhs)

    # Primal endpoint sandwich L <= Q <= P <= U.
    for w in WEIGHTS:
        add(f"Q_ge_L_{w}", {vq(w): 1}, "ge", L[w])
        add(f"P_le_U_{w}", {vq(w): 1, vh(w): 255}, "le", U[w])

    # Dual endpoint sandwich:
    #   U^perp <= P^perp <= Q^perp <= L^perp
    #
    # Let S_q(j) = sum q_w K_j(w) over full symmetric spectrum,
    #     S_h(j) = sum h_w K_j(w).
    #
    # Bq_j = S_q / 2^123
    # Bp_j = (S_q + 255 S_h) / 2^131
    # Bq_j - Bp_j = 255 * (S_q - S_h) / 2^131.
    #
    # Hence Bq >= Bp is simply S_q >= S_h.
    divisibility = []
    for j in DUAL_WEIGHTS:
        Sq = {}
        Sh = {}
        for w in WEIGHTS:
            a = symmetric_kraw_coeff(N, j, w)
            add_coeff(Sq, vq(w), a)
            add_coeff(Sh, vh(w), a)

        bp_num = dict(Sq)
        for name, a in Sh.items():
            add_coeff(bp_num, name, 255 * a)

        g_num = dict(Sq)
        for name, a in Sh.items():
            add_coeff(g_num, name, -a)

        add(f"Bp_ge_Udual_{j}", bp_num, "ge", Ud[j] * (1 << KP))
        add(f"Bq_le_Ldual_{j}", Sq, "le", Ld[j] * (1 << KQ))
        add(f"Bq_ge_Bp_{j}", g_num, "ge", 0)

        divisibility.append({
            "name": f"Bq_integral_{j}",
            "linear_form": {k: str(v) for k, v in Sq.items()},
            "modulus": str(1 << KQ),
        })
        divisibility.append({
            "name": f"Bp_integral_{j}",
            "linear_form": {k: str(v) for k, v in bp_num.items()},
            "modulus": str(1 << KP),
        })
        divisibility.append({
            "name": f"G_integral_{j}",
            "linear_form": {k: str(v) for k, v in g_num.items()},
            "modulus": str(1 << KP),
            "meaning": "G_j=(Bq_j-Bp_j)/255",
        })

    # Every fixed-weight shell of the affine-invariant endpoint codes P and Q
    # is a 2-design.  Export both the point and pair incidence divisibilities;
    # these are exact lattice constraints, not consequences of a real LP.
    for w in WEIGHTS[1:]:
        pairs = comb(w, 2)
        divisibility.extend([
            {
                "name": f"Q_shell_1_design_{w}",
                "linear_form": {vq(w): str(w)},
                "modulus": str(N),
                "meaning": "q_w * w is divisible by 256",
            },
            {
                "name": f"P_shell_1_design_{w}",
                "linear_form": {vq(w): str(w), vh(w): str(255 * w)},
                "modulus": str(N),
                "meaning": "p_w * w is divisible by 256",
            },
            {
                "name": f"Q_shell_2_design_{w}",
                "linear_form": {vq(w): str(pairs)},
                "modulus": str(comb(N, 2)),
                "meaning": "q_w * C(w,2) is divisible by C(256,2)",
            },
            {
                "name": f"P_shell_2_design_{w}",
                "linear_form": {vq(w): str(pairs), vh(w): str(255 * pairs)},
                "modulus": str(comb(N, 2)),
                "meaning": "p_w * C(w,2) is divisible by C(256,2)",
            },
        ])

    # Redundant but useful normalized first-shell lattice constraint.
    divisibility.append({
        "name": "h_38_multiple_128",
        "linear_form": {vh(38): "1"},
        "modulus": "128",
        "meaning": "p_38=255*h_38 and the P weight-38 shell is a 2-design",
    })

    # Useful shell congruences / design divisibility as metadata.
    metadata["integer_variables"] = [vq(w) for w in WEIGHTS] + [vh(w) for w in WEIGHTS]
    metadata["derived_congruences"] = [
        "p_w-q_w = 255 h_w",
        "c_w = q_w + 31 h_w",
        "At w=38: p_38=32640*m, h_38=128*m, c_38=3968*m, m>=230.",
    ]
    metadata["design_divisibility"] = {
        "description": "Each fixed-weight shell of affine-invariant P,Q is a 2-design.",
        "conditions": [
            "A_w * w is divisible by 256.",
            "A_w * C(w,2) is divisible by C(256,2)=32640.",
        ],
        "exported_for": ["Q", "P"],
    }

    out = {
        "metadata": metadata,
        "constraints": constraints,
        "divisibility_constraints": divisibility,
    }

    out_path = Path(__file__).resolve().parents[1] / "generated" / "coupled_lp_exact.json"
    out_path.parent.mkdir(exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2))
    print(f"Wrote {out_path}")
    print(f"{len(metadata['variables'])} variables")
    print(f"{len(constraints)} linear constraints")
    print(f"{len(divisibility)} divisibility constraints")


if __name__ == "__main__":
    main()
