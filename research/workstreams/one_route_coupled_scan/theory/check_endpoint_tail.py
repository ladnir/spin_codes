"""Exact bounded empty-gap tails before active events.

Retain the actual reset/live terminal function for a fixed input mask before
collapsing to the common two-state envelope. This is a local diagnostic,
not a route-averaged or whole-code distance certificate. No file output.
"""

from fractions import Fraction as F

from check_empty_span_potential import check
from check_two_groups import fresh_terms
from coupled_rank_transfer import local_matrix, pair_antipodal_cycle
from empty_span_rank import empty_step
from paired_rank_step import active_step


def active_terminals(vector, mask, z, tau, *, fresh=False, ranked=False, outer_randomized=False):
    """Two columns for one ordinary active step, on the SAME global quotient."""
    cycle = pair_antipodal_cycle(4)
    if outer_randomized:
        if fresh or ranked:
            raise ValueError("lazy outer uses its own charge, not source freshness or fixed-input ranks")
        from lazy_outer_step import terminals
        return terminals(vector, mask, z, tau)
    if ranked:
        quiet_terminal = [[F(0) for _ in range(16)] for _ in range(5)]
        quiet_terminal[0][0] = F(1)
        reset = active_step(quiet_terminal, 4, mask, z, tau, cycle=cycle, fresh=fresh)
        live = active_step(vector, 4, mask, z, tau, cycle=cycle, fresh=fresh)
        return reset, live
    matrix = local_matrix(4, mask, z, min(tau, F(256)), cycle)
    k = mask.bit_count()
    reset = [[F(0) for _ in range(16)] for _ in range(5)]
    live = [[F(0) for _ in range(16)] for _ in range(5)]
    live[0][0] = z**k * (vector[k][mask] if fresh
                         else max(vector[r][mask] for r in range(1, k + 1)))
    for h in range(1, 16):
        for r in range(1, h.bit_count() + 1):
            if fresh:
                terms = fresh_terms(mask, h, z, tau)
                reset[r][h] = sum((weight for following, weight, _ in terms if not following), F(0))
                live[r][h] = sum((weight * max(vector[s][following]
                                              for s in range(1, min(r + k - gained, following.bit_count()) + 1))
                                  for following, weight, gained in terms if following), F(0))
                continue
            reset[r][h] = matrix[h][0]
            live[r][h] = sum((matrix[h][following]
                              * max(vector[s][following]
                                    for s in range(1, min(r + k, following.bit_count()) + 1))
                              for following in range(1, 16)), F(0))
    return reset, live


def endpoint_tails(vector, mask, z, tau, steps=16, *, fresh=False, ranked=False, outer_randomized=False):
    """For each L<=steps, a linear envelope for L empty steps, then one event.

    If the final table is K_ell, a gap L>=ell is bounded by
    diag(1,lambda**(L-ell)) K_ell. Keep L<ell separately.
    """
    if not 0 <= steps <= 32:
        raise ValueError("bounded exact tail accepts 0..32 steps")
    reset, live = active_terminals(vector, mask, z, tau, fresh=fresh, ranked=ranked,
                                   outer_randomized=outer_randomized)
    cycle = pair_antipodal_cycle(4)
    valid = [(r, h) for h in range(1, 16) for r in range(1, h.bit_count() + 1)]
    tables = []
    for length in range(steps + 1):
        tables.append([[F(0), live[0][0]],
                       [max(reset[r][h] / vector[r][h] for r, h in valid),
                        max(live[r][h] / vector[r][h] for r, h in valid)]])
        if length != steps:
            reset = empty_step(reset, 4, z, tau, cycle=cycle)
            live = empty_step(live, 4, z, tau, cycle=cycle)
    return tables


def main():
    vector, z, tau, contraction = check(F(999, 1000), F(4981, 5000), tau=F(65536))
    for mask in (1, 3, 5):
        tables = endpoint_tails(vector, mask, z, tau)
        for length in (0, 1, 2, 4, 8, 16):
            ratios = [tables[length][1][j] / (tables[0][1][j] * contraction**length)
                      for j in range(2)]
            # Both comparisons are exact; decimal ratios are diagnostics only.
            print(f"mask={mask:04b} L={length}: normalized reset={float(ratios[0]):.9g}, "
                  f"live={float(ratios[1]):.9g}; "
                  f"reset_improves={ratios[0] < 1}, live_improves={ratios[1] < 1}")
    print("Exact local tails only. Short gaps must remain separate in route averages.")


if __name__ == "__main__":
    main()
