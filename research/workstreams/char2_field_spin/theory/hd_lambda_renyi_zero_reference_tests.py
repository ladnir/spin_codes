from fractions import Fraction as F
from itertools import product
from unittest.mock import patch
import unittest

from flint import arb, ctx

import hd_lambda_renyi_zero_reference as ref


def mul(a, b):
    return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2)) for i in range(2) for j in range(2))


def power(a, n):
    out = (F(1), F(0), F(0), F(1))
    for _ in range(n):
        out = mul(out, a)
    return out


def table(rows):
    return {(m, f): tuple(x if isinstance(x, arb) else ref.A.uniform.rational(x) for x in row)
            for m, row in enumerate(rows) for f in range(m//2+1)}


class ZeroReferenceTests(unittest.TestCase):
    def test_shared_column_flags_and_noncommuting_word(self):
        rows = [(F(1), F(0), F(0), F(1, 3)),
                (F(1, 3), F(1, 5), F(1, 7), F(1, 2)),
                (F(1, 9), F(1, 8), F(1, 6), F(1, 4))]
        self.assertNotEqual(mul(rows[1], rows[2]), mul(rows[2], rows[1]))
        eta, v, r, columns, sites = F(2, 5), F(3, 2), 2, 3, 2
        exact = F(0)
        # Enumerate all whole-column subsets, NOT independent regional flags.
        for flags in product(range(1 << r), repeat=columns):
            word, mass, degree = (F(1), F(0), F(0), F(1)), F(1), 0
            for flag in flags:
                m = flag.bit_count()
                mass *= eta**(r-m)*(1-eta)**m
                degree += m
                word = mul(word, power(rows[m], sites))
            exact += mass*v**degree*(word[0]+word[1])
        with ctx.workprec(256):
            prepared = ref.prepare_reference(table(rows), r=r, regions_per_column=sites, columns=columns)
            matrix = ref.column_matrix(prepared, eta, alpha=1, flag_mark=v)
            word = ref.A.uniform.scalar_power(matrix, columns)
            self.assertTrue((word[0]+word[1]).contains(ref.A.uniform.rational(exact)))
            got = ref.reference_log(prepared, eta, alpha=1, flag_mark=v, minimum_flags=4, source_count=7)
            self.assertTrue(got.contains(ref.A.log2(ref.A.uniform.rational(7*exact/v**4))))

    def test_fractional_flag_weights_not_unpowered(self):
        with ctx.workprec(256):
            prepared = ref.prepare_reference(table([(F(1), 0, 0, F(1, 4)), (F(1, 2), 0, 0, F(1, 3))]),
                                             r=1, regions_per_column=1, columns=1)
            half = ref.column_matrix(prepared, F(1, 2), alpha=F(1, 2))[0]
            one = ref.column_matrix(prepared, F(1, 2), alpha=1)[0]
            self.assertTrue(half > one)
            self.assertTrue(half.overlaps(arb(2).sqrt()*3/4))

    def test_zero_region_evolution_and_uniform_endpoint(self):
        rows = [(F(1), F(0), F(0), F(1, 8)), (F(1, 4), F(1, 3), F(1, 5), F(1, 2))]
        with ctx.workprec(256):
            prepared = ref.prepare_reference(table(rows), r=1, regions_per_column=2, columns=3)
            zero_column = ref.column_matrix(prepared, 1)
            self.assertTrue(zero_column[3].contains(arb(1)/64))
            uniform_word = ref.A.uniform.scalar_power(ref.column_matrix(prepared, 0), 3)
            for actual, expected in zip(uniform_word, power(rows[1], 6)):
                self.assertTrue(actual.contains(ref.A.uniform.rational(expected)))
            bad = table(rows)
            del bad[0, 0]
            with self.assertRaises(ValueError):
                ref.prepare_reference(bad, r=1)

    def test_positive_flag_chernoff_dominates_exact_tail(self):
        # Scalar matrices isolate the degree tail, including gapped coefficients.
        rows = [(F(1), 0, 0, F(1)), (F(1, 10), 0, 0, F(1, 10)), (F(1, 3), 0, 0, F(1, 3))]
        eta, v = F(1, 3), F(2)
        exact_tail = F(0)
        for flags in product(range(4), repeat=3):
            degree = sum(flag.bit_count() for flag in flags)
            if degree < 4:
                continue
            exact_tail += product_fraction(eta**(2-f.bit_count())*(1-eta)**f.bit_count()
                                           *rows[f.bit_count()][0] for f in flags)
        with ctx.workprec(256):
            prepared = ref.prepare_reference(table(rows), r=2, regions_per_column=1, columns=3)
            bound = ref.reference_log(prepared, eta, alpha=1, flag_mark=v, minimum_flags=4)
            self.assertTrue(arb(2)**bound > ref.A.uniform.rational(exact_tail))

    def test_exact_toy_reference_moment_atoms(self):
        # Uniform repetition code (x,x) in F4^2, then independent NZ output diagonals.
        q, width, a, eta = 4, 2, 1, F(2, 7)
        probability = {(0, 0): F(1, 4)}
        probability.update({y: F(1, 12) for y in product(range(1, 4), repeat=2)})
        for p in (2, 3, 4):
            old = sum((mass**p/F(1, 16)**(p-1) for mass in probability.values()), F(0))
            zero, nonzero = ref.reference_moment_atoms(a, p, old, eta, q=q, width=width)
            direct_zero, direct_nonzero = F(0), F(0)
            for y, mass in probability.items():
                reference = eta*(y == (0, 0))+(1-eta)*F(1, 16)
                term = mass**p/reference**(p-1)
                if y == (0, 0):
                    direct_zero += term
                else:
                    direct_nonzero += term
            self.assertEqual((zero, nonzero), (direct_zero, direct_nonzero))
            old_atoms = ref.reference_moment_atoms(a, p, old, 0, q=q, width=width)
            self.assertEqual(sum(old_atoms), old)
        full = ref.reference_moment_atoms(16, 3, 1, eta)
        self.assertNotEqual(sum(full), 1)  # P=U does NOT imply L_(Qeta)=1.

    def test_scheduled_weight_eta_zero_reproduces_frozen_activity(self):
        moments = {(): {p: F(1) for p in (2, 3, 4)}}
        with ctx.workprec(256):
            for p in (2, 3, 4):
                got, info = ref.scheduled_weights(moments, p, F(0), arities=(16,))
                old, _ = ref.activity.scheduled_weights(moments, p, arities=(16,))
                self.assertTrue(got[0].overlaps(old[0]))
                self.assertEqual(info[0]['labelled_rowsets'], 1)

    def test_exact_endpoint_max_and_saved_roundtrip(self):
        with ctx.workprec(256):
            tiny = arb(2)**-1000
            rows = table([(1, 0, 0, tiny), (tiny, tiny, tiny, tiny), (tiny, tiny, tiny, tiny)])
            rows[2, 1] = (tiny*2, tiny*3, tiny*4, tiny*5)
            prepared = ref.prepare_reference(rows, r=2)
            self.assertEqual(prepared['regional_envelope_upper'][0][1], 0)
            self.assertEqual(prepared['regional_envelope_upper'][2][0], F(1, 2**999))
            saved = [[str(x) for x in row] for row in prepared['regional_envelope_upper']]
            restored = ref.reference_from_envelope(saved, r=2)
            self.assertEqual(restored['regional_envelope_upper'], prepared['regional_envelope_upper'])

    def test_mock_receipt_provider_once_and_no_certificate_claim(self):
        fake = dict(controls=[], source_count='123', source_profiles=99)
        with patch.object(ref.A.uniform, 'regional_directed', return_value='bank') as provider, \
                patch.object(ref, 'controls_from_regions', return_value=fake) as controls:
            receipt = ref.receipt(256, flag_marks=(F(1), F(2)))
        provider.assert_called_once_with(ref.MARK, ref.ALPHA, 18, 60)
        controls.assert_called_once_with('bank', flag_marks=(F(1), F(2)))
        self.assertFalse(receipt['failure_bound'])
        self.assertFalse(receipt['whole_code_certificate'])
        self.assertTrue(receipt['whole_column_flags'])
        self.assertTrue(receipt['zero_input_regional_evolution_retained'])

    def test_source_upper_count_once_and_four_fixed_etas(self):
        profiles = {(i, 0): i+1 for i in range(99)}
        terms = [(0, 0, 1, arb(1))]
        regions = table([(1, 0, 0, F(1, 2))]*19)
        with ctx.workprec(256), patch.object(ref.A, 'profile_counts', return_value=profiles), \
                patch.object(ref.A, 'source_profile_terms', return_value=terms), \
                patch.object(ref, 'reference_log', return_value=arb(-7)) as evaluate:
            data = ref.controls_from_regions(regions)
        self.assertEqual(data['source_count'], str(sum(profiles.values())))
        self.assertEqual([row['eta'] for row in data['controls']], list(map(str, ref.ETAS)))
        self.assertEqual(evaluate.call_count, 5)  # Four eta controls plus max-f U control.
        for call in evaluate.call_args_list:
            self.assertEqual(call.kwargs['source_count'], sum(profiles.values()))

    def test_guards_and_pins(self):
        self.assertEqual(ref.dependencies()['hd_lambda_renyi_activity.py'], ref.PINS['hd_lambda_renyi_activity.py'])
        with self.assertRaises(ValueError):
            ref.reference_moment_atoms(12, 3, 1, 1)
        with self.assertRaises(ValueError):
            ref.reference_moment_atoms(1, 3, 1, F(1, 2))
        with self.assertRaises(ValueError):
            ref.scheduled_weights({}, 3, F(1, 2), arities=(10,))
        with self.assertRaises(ValueError):
            ref.receipt(128)


def product_fraction(values):
    result = F(1)
    for value in values:
        result *= value
    return result


if __name__ == '__main__':
    unittest.main()
