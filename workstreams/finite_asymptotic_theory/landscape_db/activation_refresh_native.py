"""Native four-state contraction; epoch powering uses the checked Python code."""
import ctypes
import math
from pathlib import Path

import numpy as np
import activation_q1_refresh as reference

LIBRARY = Path(__file__).resolve().parent/'native_build/activation_refresh_kernel.dll'


class RefreshKernel:
    def __init__(self):
        self.library = ctypes.CDLL(str(LIBRARY))
        self.function = self.library.refresh_log_coefficients
        pointer = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
        self.function.argtypes = [pointer,pointer,ctypes.c_int,ctypes.c_int,pointer]
        self.function.restype = ctypes.c_int

    def coefficients(self,zero,one,epochs,block):
        zero = np.ascontiguousarray(zero,dtype=np.float64)
        one = np.ascontiguousarray(one,dtype=np.float64)
        if (zero.ndim!=3 or zero.shape[1:]!=(4,4) or one.shape!=zero.shape
                or not len(zero) or not 1<=block<=4096 or epochs<1
                or any(np.isnan(m).any() or np.isposinf(m).any() for m in (zero,one))):
            raise ValueError('invalid four-state geometry or matrix')
        rz,ra = reference.region_logs(zero,one,epochs)
        ra -= math.log(epochs)
        output = np.empty((len(zero),block+1),dtype=np.float64)
        status = self.function(rz,ra,len(zero),block,output)
        if status:
            raise ArithmeticError(f'four-state recurrence failed: {status}')
        output -= np.array([math.log(math.comb(block,w)) for w in range(block+1)])
        return output
