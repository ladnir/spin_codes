"""Normalize regional coefficients before maximizing over support boxes.

For each parity, enumerate only its projected (r,f) histogram domain.
Distribute outer density and cyclic-run tilts across the 256 regions.
Positive matrix domination then bounds all profiles in the box, without
separately maximizing a coefficient numerator and its cancelling scalar.
"""
from fractions import Fraction as F

from flint import arb, ctx

from hd_shifted_alias_certificate import rational
from hd_shifted_cover_backend import validate_box, profile_count
from hd_shifted_coverage import rational_ceiling
from hd_shifted_run_count import transfer_value
from hd_shifted_unoriented_certificate import matrix_product
from hd_sumrank_paired_certificate import scalar_power


class ProjectedBound:
    def __init__(self, witness, x, y):
        # transfer_value rejects floating markers and checks positivity.
        self.count_value = transfer_value(4096, x, y)
        self.witness, self.x, self.y = witness, F(x), F(y)
        self.cache = ({}, {})
        with ctx.workprec(witness.precision):
            log2 = lambda v: v.log()/arb(2).log()
            self.log_x, self.log_y = log2(rational(self.x)), log2(rational(self.y))

    def regional(self, parity, r, f):
        """Upper-enclose R_parity(r,f) times its density/count tilt."""
        if parity not in (0,1) or any(not isinstance(v,int) for v in (r,f)):
            raise ValueError('Parity and histogram counts must be integers')
        hist = 4096-r+f, r-2*f, f
        if min(hist) < 0:
            raise ValueError('Nonnegative regional histogram required')
        key = r, f
        if key not in self.cache[parity]:
            w = self.witness
            with ctx.workprec(w.precision):
                if w.grid is None or w.uniform_aliases:
                    numerator = w.regions[parity]
                else:
                    numerator = w.regions[parity].numerator(((hist[1],hist[1]),(f,f)))
                normalizer = w.log_factorial(4096)-sum((w.log_factorial(n) for n in hist),arb(0))
                normalizer += sum((rational(n)*m for n,m in zip(hist[1:],w.log_markers[parity])),arb(0))
                # Over128 regions of each parity, these tilts multiply to
                # beta^(alpha*r) / (x^r y^(r-a-b)).
                tilt = ((w.slope-self.log_x-self.log_y)*r/256
                        +self.log_y*f/128-normalizer)
                scale = arb(2)**tilt
                self.cache[parity][key] = tuple((entry*scale).upper() for entry in numerator)
        return self.cache[parity][key]

    def envelopes(self, ranges):
        ranges = validate_box(ranges)
        if not profile_count(ranges):
            raise ValueError('Empty unsigned support box')
        result, counts = [], []
        for parity in (0,1):
            values = [arb(0) for _ in range(4)]
            count = 0
            for r in range(ranges[0][0],ranges[0][1]+1):
                lo, hi = ranges[parity+1]
                for f in range(max(lo,r-4096),min(hi,r//2)+1):
                    region = self.regional(parity,r,f)
                    values = [max(a,b) for a,b in zip(values,region)]
                    count += 1
            result.append(tuple(values))
            counts.append(count)
        return result, counts

    def receipt(self, ranges):
        ranges = validate_box(ranges)
        w = self.witness
        with ctx.workprec(w.precision):
            regions, counts = self.envelopes(ranges)
            value = scalar_power(matrix_product(*regions),128)
            # G(x,y) counts every indexed support with its weight/run tilt.
            # Taking a subset of that sum only decreases it. No profile
            # multiplier or second whole-box count belongs here.
            bound = ((rational(self.count_value)*(value[0]+value[1])).log()/arb(2).log()
                     +w.constant)
            ceiling = rational_ceiling(bound)
        return dict(ranges=ranges,integer_profile_triples=profile_count(ranges),
                    projected_histogram_counts=counts,log2_upper_rational=str(ceiling),
                    log2_upper_enclosure=str(bound),directed=True,
                    covers_all_even_odd_splits=True,whole_code_certificate=False,
                    construction=w.construction,field_order=256,forbidden_weight=734003,
                    support_count_method='projected-normalization-with-cyclic-run-tilt')

    def witness_record(self):
        return dict(self.witness.witness(),projected_normalization=True,
                    count_markers=[str(self.x),str(self.y)])


def checkpoint_receipts(precision=256):
    """Recompute the two fixed swapped-route boxes; never load saved bounds."""
    from hd_shifted_swap_certificate import SwapCachedBound
    specifications = (
        (((500,524),(188,212),(188,212)),
         '6231346319206105/9007199254740992',
         '3529422830148241/36028797018963968',
         ('21153727/526426284','62479461/686625565'), ('375/473','49/11825')),
        (((1012,1036),(388,412),(388,412)),
         '2172973163350011/4503599627370496',
         '3695833225349633/36028797018963968',
         ('68988939/780108641','179919979/892083004'), ('25/31','7/775')),
    )
    for ranges,z,alpha,row,count_markers in specifications:
        row = tuple(map(F,row))
        backend = SwapCachedBound(F(z),F(alpha),(row,row),grid=128,precision=precision)
        projected = ProjectedBound(backend,*map(F,count_markers))
        receipt = projected.receipt(ranges)
        receipt['witness'] = projected.witness_record()
        yield receipt


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision',type=int,default=256)
    args = parser.parse_args()
    for receipt in checkpoint_receipts(args.precision):
        print(json.dumps(receipt),flush=True)
