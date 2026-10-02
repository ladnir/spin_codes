"""Cross-check the research scalar oracle against the frozen proof constructors."""
import re
import unittest
from pathlib import Path

import packet_inner_small_extension as small
import packet_inner_quadratic_extension as quadratic


ROOT = Path(__file__).resolve().parent


class ImplementationMapTests(unittest.TestCase):
    def test_literal_base_and_border_match_proof(self):
        source = (ROOT/'implementation/rs16x8_border/RsBorderScalar.cpp').read_text()
        literal = re.search(r'columns\[64\]\s*=\s*\{([^}]+)\}', source).group(1)
        base = [int(v.strip(), 0) for v in literal.split(',') if v.strip()]
        self.assertEqual(len(base), 64)
        for bits in range(17, 21):
            module = small if bits <= 18 else quadratic
            _, expected, record = module.construction(bits)
            actual = [value | sum((((p&1)*((p>>j)&1)) << (15+j))
                                  for j in range(1, bits-15))
                      for p, value in enumerate(base)]
            self.assertEqual(actual, list(expected), f's{bits} implementation/proof map mismatch')
            self.assertEqual(record['appended_monomials'], [[0,j] for j in range(1,bits-15)])
        other = (ROOT/'implementation/rs16x8/RsWideScalar.cpp').read_text()
        literal = re.search(r'columns\[64\]\s*=\s*\{([^}]+)\}', other).group(1)
        self.assertEqual(base, [int(v.strip(), 0) for v in literal.split(',') if v.strip()])


if __name__ == '__main__':
    unittest.main()
