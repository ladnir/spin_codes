"""Tiny exact tests; no benchmark or whole-code distance claim."""
import unittest
from collections import Counter
from itertools import combinations, product
from random import Random

import independent_checks as check


class AlgebraTests(unittest.TestCase):
    def test_all_byte_adjoint_matrices(self):
        for scalar in range(256):
            matrix = check.gfni_adjoint_matrix(scalar)
            for bit in range(8):
                self.assertEqual(check.gfni_affine(1 << bit, matrix),
                                 check.adjoint_multiply(scalar, 1 << bit))
                for out_bit in range(8):
                    self.assertEqual(check.dot(check.multiply(scalar, 1 << bit), 1 << out_bit),
                                     check.dot(1 << bit, check.adjoint_multiply(scalar, 1 << out_bit)))
        self.assertNotEqual(check.multiply(2, 1), check.adjoint_multiply(2, 1))

    def test_local_maps(self):
        ac, cc = check.local_maps()
        self.assertEqual(check.rank(ac), 16)
        self.assertEqual(check.rank(cc), 16)
        self.assertTrue(all(check.apply_columns(cc, image) == 0 for image in ac))
        for i in range(8):
            self.assertEqual(check.rank(cc[8 * i:8 * (i + 1)]), 8)
            for j in range(i):
                self.assertEqual(check.rank(cc[8 * i:8 * (i + 1)] + cc[8 * j:8 * (j + 1)]), 16)
        # Field transpose is not binary transpose in this polynomial basis.
        self.assertNotEqual(tuple(check.apply_transpose(ac, 1 << i) for i in range(64)), cc)

    def test_gl2_transitivity_small_field_exact(self):
        mul = lambda a, b: check.multiply(a, b, degree=2, modulus=0b111)
        matrices = [m for m in product(range(4), repeat=4) if mul(m[0], m[3]) ^ mul(m[1], m[2])]
        self.assertEqual(len(matrices), (16 - 1) * (16 - 4))
        for state in range(1, 16):
            a, b = state & 3, state >> 2
            law = Counter((mul(m[0], a) ^ mul(m[1], b)) |
                          ((mul(m[2], a) ^ mul(m[3], b)) << 2) for m in matrices)
            self.assertEqual(law, Counter({target: 12 for target in range(1, 16)}))

    def test_gl2_binary_adjoint(self):
        random = Random(909)
        for _ in range(30):
            matrix = check.sample_matrix(random)
            columns = tuple(check.matrix_action(matrix, 1 << bit) for bit in range(16))
            self.assertEqual(check.rank(columns), 16)
            for bit in range(16):
                self.assertEqual(check.matrix_adjoint(matrix, 1 << bit),
                                 check.apply_transpose(columns, 1 << bit))

    def test_multistep_inner_adjoint_and_boundary(self):
        random = Random(411)
        for count in (1, 2, 3, 7, 17):
            matrices = [check.sample_matrix(random) for _ in range(count)]
            x = [random.getrandbits(64) for _ in range(count)]
            y = [random.getrandbits(64) for _ in range(count)]
            encoded, reversed_y = check.inner_forward(x, matrices), check.inner_transpose(y, matrices)
            self.assertEqual(sum(check.dot(a, b) for a, b in zip(encoded, y)) % 2,
                             sum(check.dot(a, b) for a, b in zip(x, reversed_y)) % 2)
            changed = matrices[:-1] + [check.sample_matrix(random)]
            self.assertEqual(check.inner_forward(x, matrices), check.inner_forward(x, changed))
            self.assertEqual(check.inner_transpose(y, matrices), check.inner_transpose(y, changed))

    def test_true_byte_route_and_padding(self):
        for groups in (8, 16, 24):
            route = check.byte_route(groups, Random(groups))
            self.assertEqual(len(set(route)), 32 * groups)
            expected = {260 * group + 8 * packet for group in range(groups) for packet in range(32)}
            self.assertEqual(set(route), expected)
            for region in range(32):
                self.assertEqual(sorted(base // 260 for base in route[groups * region:groups * (region + 1)]),
                                 list(range(groups)))
            self.assertTrue(all(base % 260 <= 248 for base in route))

    def test_native_outer_coordinate_cancellation(self):
        def permute(value):
            return sum(((value >> c) & 1) << (4 * (c % 4) + c // 4) for c in range(16))
        random = Random(927)
        for _ in range(8):
            scalar, coded = random.randrange(1, 65536), random.randrange(65536)
            columns = [check.tower16_multiply(scalar, 1 << c) for c in range(16)]
            rows = [sum(((columns[c] >> r) & 1) << c for c in range(16)) for r in range(16)]
            stored_reverse_rows = list(map(permute, rows))
            logical_forward = 0
            for column in range(16):
                if (coded >> column) & 1:
                    logical_forward ^= stored_reverse_rows[column]
            explicit_physical = permute(logical_forward)
            self.assertEqual(explicit_physical, check.apply_transpose(columns, coded))

    def test_packed_payload_orientation(self):
        for coordinate in range(8):
            for byte in range(16):
                for bit in range(8):
                    raw = bytearray(128)
                    raw[16 * coordinate + byte] = 1 << bit
                    packed = check.pack_payload_bytes(raw)
                    expected = bytearray(128)
                    expected[64 * (byte // 8) + 8 * (byte % 8) + 7 - bit] = 1 << coordinate
                    self.assertEqual(packed, bytes(expected))
                    self.assertEqual(check.unpack_payload_bytes(packed), bytes(raw))

    def test_half_payload_outer_store_selectors(self):
        # Compare the native two-input VBMI store with both independently
        # masked single-input stores. Distinct integer tags cover every byte.
        selectors = (
            (0x3830282018100800, 0x7870686058504840, 0x3931292119110901, 0x7971696159514941,
             0x3a322a221a120a02, 0x7a726a625a524a42, 0x3b332b231b130b03, 0x7b736b635b534b43),
            (0x3c342c241c140c04, 0x7c746c645c544c44, 0x3d352d251d150d05, 0x7d756d655d554d45,
             0x3e362e261e160e06, 0x7e766e665e564e46, 0x3f372f271f170f07, 0x7f776f675f574f47),
        )
        indices = [list(b''.join(v.to_bytes(8, 'little') for v in row)) for row in selectors]
        values = [[[64 * (8 * plane + symbol) + byte for byte in range(64)]
                   for symbol in range(8)] for plane in range(4)]
        reference, compact, written = [None] * 2048, [None] * 2048, [0] * 2048
        for plane in range(2):
            for symbol in range(8):
                pair = values[2 * plane][symbol] + values[2 * plane + 1][symbol]
                for row, selector in enumerate(indices):
                    offset = 16 * (64 * plane + 32 * row + 4 * symbol)
                    reference[offset:offset + 64] = [pair[index] for index in selector]
                    for half in range(2):
                        one = values[2 * plane + half][symbol]
                        mask = (0x55, 0xaa)[half]
                        for byte, index in enumerate(selector):
                            if mask & (1 << (byte // 8)):
                                compact[offset + byte] = one[index & 63]
                                written[offset + byte] += 1
        self.assertEqual(compact, reference)
        self.assertEqual(written, [1] * 2048)

    def test_prefetch_reverse_schedule_boundaries(self):
        for packets in (8, 16, 32):
            lookahead = packets // 8
            for count in range(1, 20):
                epoch, visited = count, []
                while epoch > lookahead:
                    epoch -= 1
                    target = 8 * epoch - packets
                    self.assertGreaterEqual(target, 0)
                    self.assertLess(target + 7, 8 * epoch)
                    self.assertLess(target + 7, 8 * count)
                    visited.append(epoch)
                while epoch:
                    epoch -= 1
                    visited.append(epoch)
                self.assertEqual(visited, list(reversed(range(count))))
            for count in range(1, 40):
                packet, visited = count, []
                while packet > packets:
                    packet -= 1
                    self.assertGreaterEqual(packet - packets, 0)
                    self.assertLess(packet - packets, packet)
                    visited.append(packet)
                while packet:
                    packet -= 1
                    visited.append(packet)
                self.assertEqual(visited, list(reversed(range(count))))

    def test_scaled_evaluation_maps(self):
        scales = (1, 3, 5, 15, 17, 51, 85, 255)
        summary = check.local_summary(scales=scales)
        self.assertTrue(summary['ca_zero'])
        self.assertEqual(summary['a_rank'], 16)
        self.assertEqual(summary['c_rank'], 16)
        self.assertEqual(summary['one_packet_ranks'], [8] * 8)
        self.assertEqual(summary['two_packet_ranks'], [16])
        self.assertEqual(min(w for w in summary['expansion_spectrum'] if w), 17)
        self.assertEqual(summary['expansion_spectrum'][17], 3)

    def test_three_packet_kernel(self):
        # For distinct field points, the two parity checks have a one-byte
        # kernel on three packets. Every nonzero kernel word has three
        # nonzero labels, hence 255/(255^3)=1/65025 cancellation probability.
        for i, j, k in combinations(range(8), 3):
            words = set()
            for value in range(1, 256):
                labels = (check.multiply(value, j ^ k),
                          check.multiply(value, i ^ k),
                          check.multiply(value, i ^ j))
                self.assertTrue(all(labels))
                word = sum(label << (8 * position) for label, position in zip(labels, (i, j, k)))
                self.assertEqual(check.feedback(word), 0)
                words.add(word)
            self.assertEqual(len(words), 255)

    def test_local_closure_lower_sandwich_tiny(self):
        # GF4, four two-bit packets, two-symbol state; exhaustive inputs.
        mul = lambda a, b: check.multiply(a, b, degree=2, modulus=0b111)
        def expand(state):
            a, b = state & 3, state >> 2
            return sum((a ^ mul(i, b)) << (2 * i) for i in range(4))
        def feedback(word):
            a = b = 0
            for i in range(4):
                x = (word >> (2 * i)) & 3
                a ^= x
                b ^= mul(i, x)
            return a | (b << 2)
        for z in (.8, .99):
            eta = 1/((1+z)**2-1)
            for state in range(1, 16):
                for occupancy in range(5):
                    total = 0.
                    syndrome_mass = [0.] * 16
                    for word in range(256):
                        if sum(bool((word >> (2 * i)) & 3) for i in range(4)) != occupancy:
                            continue
                        weight = z**(word ^ expand(state)).bit_count()
                        total += weight
                        syndrome_mass[feedback(word)] += weight
                    for target in range(16):
                        actual = (total-syndrome_mass[target])/15
                        upper = total/15 if target or occupancy else 0.
                        self.assertGreaterEqual(actual+1e-12, (1-eta)*upper)
                    if occupancy:
                        return_upper = total/15
                        self.assertGreaterEqual((total-syndrome_mass[0])/15+1e-12,
                                                (1-eta*eta)*return_upper)

    def test_whole_rs_outer_route_inner_adjoint(self):
        random = Random(613)
        groups = 8  # Each of the 32 regions is exactly one physical step.
        outer = [check.outer_columns([random.randrange(1, 65536) for _ in range(16)])
                 for _ in range(groups)]
        self.assertTrue(all(check.rank(columns) == 128 for columns in outer))
        route = check.byte_route(groups, random)
        matrices = [check.sample_matrix(random) for _ in range(4 * groups)]
        for _ in range(3):
            message, word = random.getrandbits(128 * groups), random.getrandbits(256 * groups)
            self.assertEqual(check.dot(check.code_forward(message, outer, route, matrices), word),
                             check.dot(message, check.code_transpose(word, outer, route, matrices)))


if __name__ == '__main__':
    unittest.main()
