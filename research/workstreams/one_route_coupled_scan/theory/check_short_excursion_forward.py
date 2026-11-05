#!/usr/bin/env python3
"""Exact partial certificate retaining the outer activation constraint.

SHORT_EXCURSION_FORWARD_GATE.md proves the bound. Cover all marked layouts,
all source lines, and all masks for 1<=v<=512 and at least ceil(v/4)
guaranteed closed length-two excursions. This is not a full-code proof.

The root owns runs: python -B, timeout60s,256MiB. Integer arithmetic only,
no output files, no raw large integers, no unbounded search. Each class
gets one analytically chosen row threshold; no floating optimization.
"""

from math import comb
from time import monotonic


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def ceil_div(x, y):
    return (x+y-1)//y


def main():
    started = monotonic()
    n, b, B, k, p, bits = 16384, 128, 64, 32, 4, 32
    G, region_batches = n*b//(p*B), n//B
    max_v, term_bits = 512, 80
    binomial_region = comb(B,k+1)
    group_choices = [comb(G,r) for r in range(b//p+2)]
    checked = excluded = 0
    minimum_excess_bits = 10**9
    worst = None

    for v in range(1,max_v+1):
        # m2=3c-v-2 is the guaranteed number of closed length-two
        # excursions. All other excursions are allowed without charge.
        lower_c = max(1,ceil_div(v+2+ceil_div(v,4),3))
        for c in range(lower_c,(v+1)//2+1):
            m2 = max(0,3*c-v-2)
            require(m2 >= ceil_div(v,4), "wedge coverage")
            if min(B,2*c+(v-2*c)//region_batches) <= k:
                excluded += 1
                continue
            layouts = comb(v-c,c-1)*comb(n-v+c,c)

            # Small-row probability <=q^[(k+m2)R0-m2(b+1)].
            # Pick the largest multiple of p that guarantees a per-class
            # inner contribution <2^-term_bits after all layout choices.
            required_q_charge = ceil_div(term_bits+layouts.bit_length(),bits)
            max_rows = (m2*(b+1)-required_q_charge)//(k+m2)
            R0 = max(0,min(b-p,(max_rows//p)*p))
            exponent = (k+m2)*R0-m2*(b+1)
            if R0:
                require(exponent < 0 and (layouts << term_bits) < (1 << (-bits*exponent)),
                        ("small-row term",v,c,R0,exponent))
            # R0=0 contains only the zero outer word, so its event is empty.

            r = R0//p+1
            # mu <=choose(B,k+1)(v/n)^(k+1). Here it is below1.
            mu_num = binomial_region*v**(k+1)
            mu_den = n**(k+1)
            require(mu_num < mu_den, "sparse-tail range exceeded")
            numerator = layouts*group_choices[r]*mu_num**r
            denominator = mu_den**r
            require((numerator << term_bits) < denominator,
                    ("large-row term failed",v,c,R0,r,
                     denominator.bit_length()-numerator.bit_length()-1))
            excess = denominator.bit_length()-numerator.bit_length()-1
            if excess < minimum_excess_bits:
                minimum_excess_bits, worst = excess,(v,c,R0,r,m2)
            checked += 1
        if monotonic()-started > 60:
            raise TimeoutError("60-second cap; incomplete coverage is not a certificate")

    require(checked > 0, "empty coverage")
    # Every checked class has two contributions below2^-80. Keep the
    # power-of-two union bound conservative and perform an integer check.
    union_numerator = 2*checked
    margin = term_bits-union_numerator.bit_length()
    require(margin >= 40 and (union_numerator << margin) < (1 << term_bits),
            "complete wedge union failed")
    print(dict(kind="exact_short_excursion_wedge_not_full_code",K=n*b//2,
               B=B,packet=p,batch=b,field_bits=bits,v_range=(1,max_v),
               class_condition="3*c-v-2 >= ceil(v/4), 2*c-1 <= v",
               checked_classes=checked,zero_guard_classes=excluded,
               term_bits=term_bits,union_numerator=union_numerator,
               union_denominator_power=term_bits,
               failure_margin_at_least_bits=margin,
               weakest_large_row_margin_bits=minimum_excess_bits,
               weakest_large_row_class=worst,
               seconds=round(monotonic()-started,3)))
    print("All masks, amplitudes and placements in this wedge are covered. Other classes remain open.")


if __name__ == "__main__":
    main()
