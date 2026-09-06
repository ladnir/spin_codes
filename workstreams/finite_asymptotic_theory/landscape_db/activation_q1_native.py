"""Native implementation of the existing log-domain Q1 recurrence."""
import ctypes
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LIBRARY = HERE / 'native_build' / 'activation_q1_kernel.dll'


class Q1Kernel:
    def __init__(self):
        self.library = ctypes.CDLL(str(LIBRARY))
        self.function = self.library.q1_log_coefficients
        pointer = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
        self.function.argtypes = [pointer, pointer, ctypes.c_int, ctypes.c_int, pointer]
        self.function.restype = ctypes.c_int

    def coefficients(self, zero, one, length):
        zero = np.ascontiguousarray(zero, dtype=np.float64)
        one = np.ascontiguousarray(one, dtype=np.float64)
        if (zero.ndim != 3 or zero.shape[1:] != (3,3) or one.shape != zero.shape
                or not len(zero) or not 1 <= length <= 4096
                or any(np.isnan(m).any() or np.isposinf(m).any() for m in (zero,one))):
            raise ValueError('invalid Q1 coefficient geometry or matrix')
        output = np.empty((len(zero), length+1), dtype=np.float64)
        status = self.function(zero, one, len(zero), length, output)
        if status:
            raise ArithmeticError(f'native Q1 recurrence failed: {status}')
        output -= np.array([math.log(math.comb(length,w)) for w in range(length+1)])
        return output
