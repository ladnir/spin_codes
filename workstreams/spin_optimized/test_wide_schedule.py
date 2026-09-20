"""Small, upstream-independent regression tests for BCH scheduling checks."""
import unittest
from wide_schedule import reorder, xor_inputs


def fixture():
    lines = ["static void bchForward(const Vec* a,Vec* x) {"]
    lines += [f"const auto o{i}=Ops::load(a+{i});" for i in range(128)]
    lines += ["const auto g0=vx(o1,o1);", "const auto z0=vx(o0,g0);"]
    lines += [f"const auto z{i}=o{i % 128};" for i in range(1, 256)]
    lines += [f"Ops::store(x+{i},z{i});" for i in range(256)]
    words = [sum(1 << (j % 64) for j in range(64*w, 64*w+64) if j % 128 == i)
             for i in range(128) for w in range(4)]
    return "\n".join(lines) + "\n}\n", ",".join(f"0x{w:x}ULL" for w in words)


class WideScheduleTests(unittest.TestCase):
    def setUp(self):
        self.source, self.matrix = fixture()

    def test_exact_and_idempotent(self):
        output, count = reorder(self.source, self.matrix)
        self.assertEqual(count, 2)
        self.assertEqual(reorder(output, self.matrix), (output, count))
        self.assertLess(output.index("Ops::store(x+0,z0);"), output.index("const auto o2="))

    def test_rejects_changed_output(self):
        with self.assertRaisesRegex(ValueError, "differs from matrix"):
            reorder(self.source.replace("z1=o1", "z1=o2"), self.matrix)

    def test_rejects_changed_matrix(self):
        with self.assertRaisesRegex(ValueError, "differs from matrix"):
            reorder(self.source, self.matrix.replace("0x1ULL", "0x0ULL", 1))

    def test_rejects_malformed_circuit(self):
        mutations = (
            self.source.replace("Ops::store(x+1,z1);", ""),
            self.source.replace("Ops::store(x+1,z1);", "Ops::store(x+0,z0);"),
            self.source.replace("Ops::store(x+1,z1);", "Ops::store(x+1,z2);"),
            self.source.replace("z1=o1", "z1=z1"),
            self.source.replace("z1=o1", "z1=g999"),
            self.source.replace("z1=o1", "z1=vx(o1)"),
            self.source.replace("z1=o1", "z0=o1"),
            self.source.replace("Ops::load(a+0)", "Ops::load(a+1)"),
            self.source.replace("\n}", "\nOps::store(x+0,o1);\n}"),
        )
        for source in mutations:
            with self.subTest(source=source[-150:]), self.assertRaises(ValueError):
                reorder(source, self.matrix)

    def test_binary_grammar(self):
        self.assertEqual(xor_inputs("vx(o0,vx(g2,z1))"), ["o0", "g2", "z1"])
        for expression in ("", "vx(o0)", "vx(o0,o1,o2)", "o0o1", "vx(o0,o1);bad()"):
            with self.subTest(expression=expression), self.assertRaises(ValueError):
                xor_inputs(expression)


if __name__ == "__main__":
    unittest.main()
