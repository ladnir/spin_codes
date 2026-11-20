"""Independent algebra checks of the fixed S5 BCH correction basis."""
import random
import unittest
import parity_corrector as corrector


def echelon(words):
    rows = {}
    for word in words:
        while word:
            pivot = word.bit_length()-1
            if pivot in rows:
                word ^= rows[pivot]
            else:
                rows[pivot] = word
                break
    return rows


def remainder(word, rows):
    for pivot in sorted(rows, reverse=True):
        if word >> pivot & 1:
            word ^= rows[pivot]
    return word


class CorrectorBasisTests(unittest.TestCase):
    def test_fixed_code_and_right_inverse(self):
        basis = corrector.fixed_c_basis()
        rows = echelon(basis)
        self.assertEqual(len(rows), 128)
        self.assertEqual(remainder((1 << 256)-1, rows), 0)
        generator = corrector.generator_polynomial(37)
        p_basis = []
        for i in range(131):
            value = generator << i
            p_basis.append(value | ((value.bit_count() % 2) << 255))
        p_rows = echelon(p_basis)
        self.assertEqual(len(p_rows), 131)
        powers = [corrector.gf_pow(corrector.gf_pow(corrector.ALPHA, i), 37)
                  for i in range(255)] + [0]
        for word in basis:
            self.assertEqual(word.bit_count() % 2, 0)
            self.assertEqual(remainder(word, p_rows), 0)
            syndrome = 0
            for i, power in enumerate(powers):
                if word >> i & 1:
                    syndrome ^= power
            self.assertLess(syndrome, 32)
        rng = random.Random(83742109)
        permutations = [rng.sample(range(256), 256) for _ in range(3)]
        rank, section = corrector.right_inverse(basis, permutations)
        self.assertEqual(rank, 255)
        for target in range(255):
            combined = 0
            for j, pi in enumerate(permutations):
                codeword = 0
                for i, word in enumerate(basis):
                    if section[target] >> (128*j+i) & 1:
                        codeword ^= word
                self.assertEqual(remainder(codeword, rows), 0)
                # Direct coordinate loop independent of corrector.permute.
                for i in range(256):
                    if codeword >> i & 1:
                        combined ^= 1 << pi[i]
            self.assertEqual(combined, (1 << target) | (1 << 255))
        # Identical reserved permutations cannot span all even vectors.
        rank, section = corrector.right_inverse(basis, [list(range(256))]*3)
        self.assertEqual(rank, 128)
        self.assertIsNone(section)


if __name__ == '__main__':
    unittest.main()
