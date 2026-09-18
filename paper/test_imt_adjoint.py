"""Exact binary check of the manuscript's column-oriented IMT recurrences.

Compare the reverse recurrence with an independently assembled full matrix.
This checks the algebra, not the native implementation or distance certificate.
"""

from pathlib import Path
import random
import re
import unittest


def apply_columns(columns, vector):
    output = 0
    while vector:
        bit = vector & -vector
        output ^= columns[bit.bit_length() - 1]
        vector ^= bit
    return output


def transpose_apply(columns, vector):
    return sum(((column & vector).bit_count() % 2) << i
               for i, column in enumerate(columns))


def forward(vector, a, c, schedule):
    width = len(c)
    mask = (1 << width) - 1
    state = output = 0
    for i, (u, v) in enumerate(schedule):
        x = (vector >> (i * width)) & mask
        y = x ^ apply_columns(a, state)
        output |= y << (i * width)
        state ^= u if (v & state).bit_count() % 2 else 0
        state ^= apply_columns(c, x)
    return output


def reverse(vector, a, c, schedule):
    width = len(c)
    mask = (1 << width) - 1
    state = output = 0
    for i in reversed(range(len(schedule))):
        u, v = schedule[i]
        x = (vector >> (i * width)) & mask
        output |= (x ^ transpose_apply(c, state)) << (i * width)
        state ^= v if (u & state).bit_count() % 2 else 0
        state ^= transpose_apply(a, x)
    return output


class IMTAdjointTests(unittest.TestCase):
    def check_maps(self, a, c, steps):
        rng = random.Random(0x494D54 + steps)
        schedule = []
        for _ in range(steps):
            u = rng.randrange(1, 1 << len(a))
            v = rng.randrange(1 << len(a))
            if (u & v).bit_count() % 2:
                v ^= u & -u
            self.assertEqual((u & v).bit_count() % 2, 0)
            schedule.append((u, v))
        size = steps * len(c)
        # Each basis output is a column of the complete forward linear map.
        full_columns = [forward(1 << i, a, c, schedule) for i in range(size)]
        for i in range(size):
            self.assertEqual(reverse(1 << i, a, c, schedule),
                             transpose_apply(full_columns, 1 << i))

    def test_independent_small_maps(self):
        self.check_maps([0b101, 0b011], [0b01, 0b11, 0b10], 4)

    def test_selected_manuscript_maps(self):
        text = (Path(__file__).parent / 'structured_imt_appendix.tex').read_text()
        rows = re.findall(
            r'^(\d+) & \\texttt\{([0-9a-f]+)\} & \\texttt\{([0-9a-f]+)\}',
            text, re.M)
        self.assertEqual([int(i) for i, _, _ in rows], list(range(19)))
        a = [int(word, 16) for _, word, _ in rows]
        ct = [int(word, 16) for _, _, word in rows]
        c = [transpose_apply(ct, 1 << i) for i in range(128)]
        for steps in (1, 3):
            with self.subTest(steps=steps):
                self.check_maps(a, c, steps)


if __name__ == '__main__':
    unittest.main()
