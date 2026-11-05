r"""Directed35% bound for actual source counts7584 through8192.

This is a class certificate for the unchanged cheap-lane construction,
not a whole-code certificate. It concerns arbitrary GF128 messages,
including full GF256 lane span, and the actual65280-element refresh set.
The frozen older dense-band helper is neither imported nor modified.

Improved input-restricted kernel bound, pointwise in route and D1:
Let q=256, g=z+(1-z)/q, I the allowed packet inputs, ell=|I|,
H the incoming state mask, h=|H|, and J the outgoing zero coordinates.
Set R=H16 D1, r_f=rank R[J,H]=min(|J|,h), and
b=rank[R[J,H], R[J,:] C_I]=r_f+d. Consider independent uniform
x in E^I and u in E^H. Their feedback constraint has probability q^-b.
On its kernel V, output Y=x+A_Hu is uniform on a linear image W.
If Y=0, then x=A_Hu and CA=0 gives R[J,H]u=0. Therefore
dim ker(Y|V)<=h-r_f, so dim W>=ell+h-b-(h-r_f)=ell-d.
A systematic projection of W is uniform on dim W coordinates. Dropping
the remaining output coordinates only increases z^wt for0<z<1. Hence

  K_I=E[z^wt(Y) 1[feedback=0]]
     <=q^-b g^(ell-d)=q^-r_f g^ell(qg)^-d<=q^-r_f g^ell.

The last inequality uses qg>=1. This does NOT assert a comparison with
the exact full-input kernel. The final packet has no feedback condition;
the same systematic projection gives its bound g^ell without a refresh.

Use the preserved rank-potential reduction with P_q(D)=(q^D-1)/(q-1).
Its fresh-pivot contraction1/q is valid because65280>=q. Fix an actual
source support T of size r. The E-defined initial kernel has a nonzero
projection onto every group in T, or that actual support is impossible.
A uniform nonzero scalar kernel vector misses each group with probability
less than1/q. The probability of at most32 missing groups is at least
c=1-r/(33q)>=1/33. This bounds the initial potential by the sum of
scalar witnesses with supports S subset T, |T\S|<=32, divided by255c.

Apply the preserved scalar outer density beta=q^160/255^32 to each
group of S. The positive output-zero identity is used before this step;
no scalar shadow is assumed to preserve the actual state support.
Every route has sum_t|I_t(S)|=256|S|. Thus the new packet bound gives
the complete source factor A^|S|, A=beta*g^256. Missing j groups costs
A^-j, not q^(256j)/beta^j. Uniformly in T, the missing-group factor is

  C32=sum_(j=0)^32 binom(8192,j) A^-j.

After counting destination masks the continuation matrix is
B[h,j]=binom(16,j)q^-min(16-j,h). For v[h]=q^-h+q^-8 and
rho=(1+1/q)^16+65536q^-8, the elementary inequality
q^-min(a,h)<=q^-a+q^-h gives Bv<=rho*v. Also v[0]/min(v)=q^8.
There are32767 charged transitions, followed by the free final packet.

Fix N=2097152, L=floor(7N/20)=734003, r0=7584, z=3/1250.
The checker verifies0<A<1. Union over at most2^8192 actual supports
with r>=r0 gives

  U=2^8192 *33/255 *A^r0*C32*q^8*rho^32767*z^-L <2^-768.

All statements are uniform in the route. The probability averages the
ideal independent setup; it is not a certificate for one benchmark seed.
Arb evaluates this fixed expression without optimization. No result file,
build, benchmark, or bulk data is created. Require20 GiB free.
"""

import argparse
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import shutil
import time


def verify(precision):
    if precision not in (192, 256):
        raise ValueError("precision must be192 or256 bits")
    started = time.monotonic()
    source = Path(__file__).resolve()
    free_bytes = shutil.disk_usage(source.parent).free
    if free_bytes < 20*(1 << 30):
        raise RuntimeError("at least20 GiB free required")
    from flint import arb, ctx
    ctx.prec = precision
    q = arb(256)
    groups, r0, packets, records = 8192, 7584, 32768, 1 << 21
    limit = 7*records//20
    beta = q**160/arb(255)**32
    z = arb(3)/1250
    g = z+(1-z)/q
    group_factor = beta*g**256
    if not bool(group_factor > 0 and group_factor < 1):
        raise ArithmeticError("required group-factor monotonicity not verified")
    missing_factor = sum(arb(comb(groups, j))*group_factor**(-j)
                         for j in range(33))
    rho = (1+1/q)**16+65536*q**(-8)
    bound = (arb(2)**groups*33/255*group_factor**r0*missing_factor
             *q**8*rho**(packets-1)*z**(-limit))
    if time.monotonic()-started > 50:
        raise RuntimeError("bounded50-second replay exceeded")
    passed = bool(bound.is_finite() and bound < arb(2)**(-768))
    if not passed:
        raise ArithmeticError("directed bound did not prove768 bits")
    return {"schema": "cheap-lane-dense-gpadding-v1",
            "scope": "actual active source counts7584 through8192",
            "whole_code_certificate": False, "full_field_messages": True,
            "directed": True, "precision_bits": precision,
            "distance": "7/20", "output_weight_limit": limit,
            "z": "3/1250", "source_groups": groups,
            "minimum_active_source_groups": r0,
            "maximum_missing_shadow_groups": 32,
            "charged_refreshes": packets-1, "target_margin_bits": 768,
            "proved_below_target": passed,
            "log2_upper_interval": str(bound.log()/arb(2).log()),
            "source_sha256": sha256(source.read_bytes()).hexdigest(),
            "free_bytes_at_start": free_bytes, "generated_file_bytes": 0,
            "elapsed_seconds": time.monotonic()-started}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--precision", type=int, choices=(192, 256), default=192)
    args = parser.parse_args()
    print(json.dumps(verify(args.precision), sort_keys=True), flush=True)
