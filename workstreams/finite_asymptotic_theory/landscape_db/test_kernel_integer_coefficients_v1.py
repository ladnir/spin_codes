"""Check the exact recurrence against independent polynomial products."""
import unittest

import kernel_integer_coefficients_v1 as exact


class KernelCoefficients(unittest.TestCase):
    def test_direct_products(self):
        for kernel in ([1, 0, 3, 0, 1], [1, 0, 0, 0, 9, 0, 2], [1, 0, 7]):
            base = kernel[::2]
            for epochs in (1, 2, 7, 32):
                expected = [1]
                for _ in range(epochs):
                    out = [0]*(len(expected)+len(base)-1)
                    for i, a in enumerate(expected):
                        for j, b in enumerate(base):
                            out[i+j] += a*b
                    expected = out
                for maximum in (0, 5, len(expected)*2+8):
                    self.assertEqual(exact.coefficients(kernel, epochs, maximum), expected[:maximum//2+1])

    def test_rejects_odd_kernel(self):
        with self.assertRaises(ValueError):
            exact.coefficients([1, 1, 1], 2, 4)


if __name__ == '__main__':
    unittest.main()
