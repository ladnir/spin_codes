"""High-precision direct DFT comparisons for the explicit outward FFT."""
import itertools
import unittest
import numpy as np
from flint import acb, arb, ctx
import outward_radix2_fft as fft


class FFTTests(unittest.TestCase):
    def test_encloses_direct_dft(self):
        ctx.prec = 512
        rng = np.random.default_rng(734902)
        for shape, initial_error in itertools.product([(8,), (4, 8), (4, 4, 4)], [0., 2**-40]):
            values = rng.normal(size=shape)+1j*rng.normal(size=shape)
            errors = np.full(shape, initial_error)
            output, bounds = fft.transform(values, errors)
            for target in np.ndindex(shape):
                exact = acb(0)
                for source in np.ndindex(shape):
                    angle = -2*arb.pi()*sum(arb(a*b)/n for a, b, n in zip(source, target, shape))
                    # Perturb each exact input within its supplied radius.
                    v = values[source]
                    original = acb(arb(float(v.real))+arb(initial_error)/4, arb(float(v.imag)))
                    exact += original*acb(0, angle).exp()
                exact /= values.size
                computed = output[target]
                delta = abs(exact-acb(float(computed.real), float(computed.imag)))
                self.assertLessEqual(delta.upper(), arb(float(bounds[target])))

    def test_square_enclosure(self):
        ctx.prec = 512
        values = np.array([0j, 1+2j, 1e-150+2e-150j, -0.01+0.002j])
        radii = np.full(len(values), 2**-45)
        for _ in range(6):
            old, old_error = values, radii
            values, radii = fft.square(values, radii)
            for i, z in enumerate(old):
                exact = (acb(float(z.real), float(z.imag))+arb(float(old_error[i]))/2)**2
                computed = acb(float(values[i].real), float(values[i].imag))
                self.assertLessEqual(abs(exact-computed).upper(), arb(float(radii[i])))


if __name__ == '__main__':
    unittest.main()
