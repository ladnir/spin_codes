"""Explicit radix-two FFT with absolute complex error at every output.

Twiddles are Arb-certified with complex error <=4u. The butterfly bound is
ea+eb+16u*|b|+2u*(|a|+|computed_twiddle_product|)+2^-1000, u=2^-53.
Only basic IEEE binary64 operations and square root are assumed.
"""
import numpy as np
from flint import arb


def up(x):
    return np.nextafter(x, np.inf)


def magnitude(z):
    return up(np.sqrt(up(up(z.real*z.real)+up(z.imag*z.imag))))


def twiddles(n, sign=-1):
    out = []
    for j in range(n//2):
        angle = sign*2*arb.pi()*j/n
        re, im = angle.cos(), angle.sin()
        x, y = float(re), float(im)
        assert abs(arb(x)-re).upper() <= arb(2)**-52
        assert abs(arb(y)-im).upper() <= arb(2)**-52
        out.append(complex(x, y))
    return np.array(out, dtype=np.complex128)


def square(z, error):
    mag = magnitude(z)
    new_error = up(up(up(up(2*mag*error)+up(error*error))+
                      up(up(mag*mag)*2**-50))+2**-1000)
    return z*z, new_error


def axis_transform(z, error, axis):
    n = z.shape[axis]
    assert n >= 2 and n & (n-1) == 0
    shifted = np.moveaxis(z, axis, -1)
    shape = shifted.shape
    work = np.array(shifted, order='C', copy=True).reshape(-1, n)
    radii = np.array(np.moveaxis(error, axis, -1), order='C', copy=True).reshape(-1, n)
    bits = n.bit_length()-1
    reverse = np.array([int(format(i, f'0{bits}b')[::-1], 2) for i in range(n)])
    roots = {length:twiddles(length) for length in [1 << k for k in range(1, bits+1)]}
    # Batch the butterflies to bound scratch memory independently of table size.
    for start in range(0, len(work), 512):
        stop = min(start+512, len(work))
        block = work[start:stop, reverse].copy()
        bounds = radii[start:stop, reverse].copy()
        for stage in range(1, bits+1):
            width = 1 << stage
            half = width//2
            view = block.reshape(-1, n//width, width)
            ev = bounds.reshape(-1, n//width, width)
            left = view[..., :half].copy()
            right = view[..., half:]
            product = right*roots[width]
            e = up(up(up(ev[..., :half]+ev[..., half:])+
                      up(magnitude(right)*2**-49))+
                   up(up(magnitude(left)+magnitude(product))*2**-52))
            e = up(e+2**-1000)
            view[..., :half] = left+product
            view[..., half:] = left-product
            ev[..., :half] = e
            ev[..., half:] = e
        work[start:stop] = block
        radii[start:stop] = bounds
    return np.moveaxis(work.reshape(shape), -1, axis), np.moveaxis(radii.reshape(shape), -1, axis)


def transform(z, error):
    assert z.shape == error.shape and np.all(error >= 0)
    for axis in range(z.ndim):
        z, error = axis_transform(z, error, axis)
        print('Outward FFT axis', axis, 'complete', flush=True)
    scale = z.size
    assert scale & (scale-1) == 0
    return z/scale, up(error/scale)
