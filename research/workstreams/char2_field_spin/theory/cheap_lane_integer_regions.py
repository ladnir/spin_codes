"""Positive integer regional arithmetic with directed dyadic rounding.

IntegerArithmetic accepts exact nonnegative rational coefficient caps.
It represents a matrix polynomial by nonnegative integer coefficients
times one shared power of two. Every product is an exact fmpz_poly matrix
product, followed by separate axis folding/truncation. Coefficients are
then rounded UP to a fixed mantissa budget, with a shared updated exponent.
Thus the representation remains a coefficientwise upper bound by positivity.
No FLOAT transform, statistical error estimate, or interval radius is used.

full and single are maximum retained indices. A cyclic axis has modulus
maximum+1; an exact axis is a downward truncation. Packing uses stride
2*full+1, so the full degree cannot carry during an ordinary product.
Keep the regular and terminal matrices' shared exponents separate.
Six updates T=R*T followed by R=R*R produce H^64 and H^63*J.

The module does not alter local caps, introduce alpha, normalize a route
coefficient, or apply a source factor/count. Rounding may enlarge a bound;
the caller must check the resulting bound rather than assume a loss budget.
For cyclic products rounded DC dominates the original DC. Exact/mixed DC
has no such one-sided comparison after truncation and rounding together.

Import is inert. Standalone execution runs exact small tests only, with
300s/512MiB/20GiB guards and stdout output. No coefficient files are written.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
MEMORY=512*(1<<20)


def dyadic(exponent):
    return F(1<<exponent) if exponent>=0 else F(1,1<<(-exponent))


def floor_log2(value):
    if value<=0: raise ValueError("positive rational required")
    n,d=value.numerator,value.denominator
    exponent=n.bit_length()-d.bit_length()
    below=(n < d<<exponent) if exponent>=0 else (n<<(-exponent) < d)
    if below:
        exponent-=1
    return exponent


def ceil_scaled(value,exponent):
    """ceil(value/2^exponent), with exact integer arithmetic."""
    if value<0: raise ValueError("nonnegative cap required")
    n,d=value.numerator,value.denominator
    if exponent>=0: d<<=exponent
    else: n<<=-exponent
    q,r=divmod(n,d)
    return q+bool(r)


class DyadicMatrix:
    __slots__=("entries","exponent","layout")

    def __init__(self,entries,exponent,layout):
        self.entries=tuple(tuple(row) for row in entries)
        self.exponent=int(exponent)
        self.layout=layout

    @property
    def shape(self):
        return len(self.entries),len(self.entries[0])


class IntegerArithmetic:
    """Exact positive convolution plus upward coefficient quantization.

    Public pack input: {(f,s):(H00,H01,H10,H11,J0,J1)} of Fractions.
    pack returns two DyadicMatrix objects. multiply accepts compatible
    matrix shapes. extract and mass return tuples of exact Fractions.
    The caller owns guard(), source/count factors and regional normalization.
    Production mantissas96/128/192 are supported; small tests use8 bits.
    """

    def __init__(self,full,single,cyclic_full=False,cyclic_single=False,bits=128,guard=lambda:None):
        if type(full) is not int or type(single) is not int or not 0<=full<=1024 or not 0<=single<=1024:
            raise ValueError("axis maxima must be integers0..1024")
        if type(bits) is not int or not 8<=bits<=256:
            raise ValueError("mantissa budget must be8..256 bits")
        self.full,self.single,self.bits=full,single,bits
        self.cyclic_full,self.cyclic_single=bool(cyclic_full),bool(cyclic_single)
        self.stride=2*full+1; self.length=single*self.stride+full+1
        if self.length>300000: raise ValueError("packed input limit300000 coefficients per entry")
        self.product_length=2*full+self.stride*(2*single)+1
        self.layout=(full,single,self.cyclic_full,self.cyclic_single)
        self.guard=guard
        from flint import fmpz_poly
        self.poly=fmpz_poly
        self.products=self.polynomial_products=self.quantizations=0
        self.initial_quantizations=self.rounded_coefficients=self.max_shift=0
        self.max_integer_bits=self.zero_objects=0
        self.max_relative_mass_increase=F(0)

    def _check(self,matrix):
        if not isinstance(matrix,DyadicMatrix) or matrix.layout!=self.layout:
            raise ValueError("matrix layout mismatch")
        if not matrix.entries or not matrix.entries[0]: raise ValueError("empty matrix shape")
        width=len(matrix.entries[0])
        if any(len(row)!=width for row in matrix.entries): raise ValueError("ragged matrix")

    def _reduce(self,f,s):
        if self.cyclic_full: f%=self.full+1
        elif f>self.full: return None
        if self.cyclic_single: s%=self.single+1
        elif s>self.single: return None
        return f+self.stride*s

    def _finish(self,arrays,exponent,maximum_bits):
        """Arrays contain only exact nonnegative product coefficients."""
        self.guard(); shift=max(0,maximum_bits-self.bits)
        if shift:
            # A ceiling can carry into bit(bits+1); reserve another bit only
            # when some retained coefficient would actually cause that carry.
            threshold=((1<<self.bits)-1)<<shift
            if any(max(coefficients,default=0)>threshold for row in arrays for coefficients in row):
                shift+=1
        self.max_shift=max(self.max_shift,shift)
        self.max_integer_bits=max(self.max_integer_bits,maximum_bits)
        self.quantizations+=1
        entries=[]
        for row in arrays:
            output=[]
            for coefficients in row:
                self.guard()
                if shift:
                    old_mass=sum(coefficients); new_mass=0
                    for i,x in enumerate(coefficients):
                        if i%4096==0: self.guard()
                        if x:
                            y=((x-1)>>shift)+1
                            self.rounded_coefficients+=(y<<shift)!=x
                            coefficients[i]=y; new_mass+=y
                    if old_mass:
                        increase=F((new_mass<<shift)-old_mass,old_mass)
                        assert increase>=0
                        self.max_relative_mass_increase=max(self.max_relative_mass_increase,increase)
                output.append(self.poly(coefficients))
            entries.append(tuple(output))
        if maximum_bits==0:
            exponent=0; self.zero_objects+=1
        else:
            exponent+=shift
        return DyadicMatrix(entries,exponent,self.layout)

    def pack_matrix(self,matrix):
        """Pack a matrix of sparse {(f,s):Fraction} coefficient dictionaries."""
        self.guard()
        if not matrix or not matrix[0] or any(len(row)!=len(matrix[0]) for row in matrix):
            raise ValueError("nonempty rectangular matrix required")
        reduced=[]; maximum=F(0)
        for row in matrix:
            output=[]
            for item in row:
                coefficients=defaultdict(F)
                for (f,s),value in item.items():
                    if type(f) is not int or type(s) is not int or min(f,s)<0:
                        raise ValueError("nonnegative integer degrees required")
                    value=F(value)
                    if value<0: raise ValueError("nonnegative rational coefficient required")
                    index=self._reduce(f,s)
                    if index is not None and value: coefficients[index]+=value
                if coefficients: maximum=max(maximum,max(coefficients.values()))
                output.append(coefficients)
            reduced.append(output)
        exponent=floor_log2(maximum)-self.bits+1 if maximum else 0
        arrays=[]; maximum_bits=0
        for row in reduced:
            output=[]
            for item in row:
                self.guard(); coefficients=[0]*self.length
                old_mass=sum(item.values(),F(0)); new_mass=0
                for index,value in item.items():
                    n=ceil_scaled(value,exponent); coefficients[index]=n; new_mass+=n
                    maximum_bits=max(maximum_bits,n.bit_length())
                    self.rounded_coefficients+=F(n)*dyadic(exponent)!=value
                if old_mass:
                    increase=F(new_mass)*dyadic(exponent)/old_mass-1
                    assert increase>=0
                    self.max_relative_mass_increase=max(self.max_relative_mass_increase,increase)
                output.append(coefficients)
            arrays.append(output)
        self.initial_quantizations+=1
        return self._finish(arrays,exponent,maximum_bits)

    def pack(self,degrees):
        """Six rational channels, with marker powers already included."""
        maps=[{} for _ in range(6)]
        for degree,row in degrees.items():
            if len(row)!=6: raise ValueError("six H/J channels required")
            for j,value in enumerate(row): maps[j][degree]=F(value)
        return (self.pack_matrix(((maps[0],maps[1]),(maps[2],maps[3]))),
                self.pack_matrix(((maps[4],),(maps[5],))))

    def multiply(self,left,right):
        """Exact matrix convolution, axis reduction, then positive ceiling."""
        self._check(left); self._check(right); self.guard()
        rows,inner=left.shape; other,cols=right.shape
        if inner!=other: raise ValueError("matrix product dimension mismatch")
        arrays=[]; maximum_bits=0; stride=self.stride
        for i in range(rows):
            output=[]
            for j in range(cols):
                self.guard(); raw=self.poly([])
                for k in range(inner):
                    raw+=left.entries[i][k]*right.entries[k][j]
                    self.polynomial_products+=1
                    self.guard()
                if len(raw)>self.product_length: raise ArithmeticError("packed carry bound exceeded")
                coefficients=[0]*self.length
                for s in range(self.single+1):
                    if s%8==0: self.guard()
                    second_s=s+self.single+1
                    fold_s=self.cyclic_single and second_s<=2*self.single
                    for f in range(self.full+1):
                        index=f+stride*s; value=raw[index]
                        second_f=f+self.full+1
                        fold_f=self.cyclic_full and second_f<=2*self.full
                        if fold_f: value+=raw[second_f+stride*s]
                        if fold_s:
                            value+=raw[f+stride*second_s]
                            if fold_f: value+=raw[second_f+stride*second_s]
                        n=int(value)
                        if n<0: raise ArithmeticError("negative exact integer product coefficient")
                        coefficients[index]=n; maximum_bits=max(maximum_bits,n.bit_length())
                del raw
                output.append(coefficients)
            arrays.append(output)
        self.products+=1
        return self._finish(arrays,left.exponent+right.exponent,maximum_bits)

    def integer_coefficient(self,matrix,row,column,f,s):
        self._check(matrix)
        if min(f,s)<0: raise ValueError("negative degree")
        index=self._reduce(f,s)
        if index is None: raise ValueError("target outside exact retained domain")
        return int(matrix.entries[row][column][index])

    def extract(self,matrix,f,s):
        self._check(matrix); scale=dyadic(matrix.exponent)
        return tuple(tuple(self.integer_coefficient(matrix,i,j,f,s)*scale
                     for j in range(matrix.shape[1])) for i in range(matrix.shape[0]))

    def mass(self,matrix):
        """Exact represented DC, including all upward rounding."""
        self._check(matrix); scale=dyadic(matrix.exponent); out=[]
        for row in matrix.entries:
            result=[]
            for poly in row:
                self.guard()
                result.append(sum((int(poly[f+self.stride*s]) for s in range(self.single+1)
                                   for f in range(self.full+1)),0)*scale)
            out.append(tuple(result))
        return tuple(out)

    def report(self):
        return {"arithmetic":"nonnegative fmpz_poly; exact fold/truncate; upward dyadic ceiling",
            "full_limit":self.full,"single_limit":self.single,
            "cyclic_full":self.cyclic_full,"cyclic_single":self.cyclic_single,
            "stride":self.stride,"mantissa_bits":self.bits,
            "input_packed_length":self.length,"maximum_product_length":self.product_length,
            "matrix_products":self.products,"polynomial_products":self.polynomial_products,
            "initial_quantizations":self.initial_quantizations,"total_quantizations":self.quantizations,
            "rounded_coefficients":self.rounded_coefficients,"maximum_shift":self.max_shift,
            "maximum_prerounding_integer_bits":self.max_integer_bits,"zero_objects":self.zero_objects,
            "maximum_entry_relative_mass_increase_per_rounding":str(self.max_relative_mass_increase),
            "error_budget_assumed":False,"alpha_applied":False,"source_factor_applied":False,
            "support_count_applied":False,"route_denominator_applied":False,"generated_bytes":0}


def exact_tests(guard=lambda:None):
    """Independent Fraction ring products, six doublings and carry probes."""
    checks=0
    carry=IntegerArithmetic(0,0,bits=8,guard=guard)
    rounded=carry._finish([[[65535]]],0,16)
    assert rounded.exponent==9 and int(rounded.entries[0][0][0])==128
    assert carry.extract(rounded,0,0)[0][0]==65536
    checks+=3
    for value in (F(0),F(1,3),F(7,16),F(3),F(999,7)):
        for e in (-9,-2,0,3,10):
            n=ceil_scaled(value,e)
            assert n*dyadic(e)>=value
            if n: assert (n-1)*dyadic(e)<value
            checks+=1
        if value:
            e=floor_log2(value); assert dyadic(e)<=value<dyadic(e+1); checks+=1
    a=(( {(0,0):F(1,3),(1,0):F(2,7)}, {(0,1):F(3,11)} ),
       ( {(1,1):F(5,13)}, {(0,0):F(7,17),(2,0):F(1,19)} ))
    b=(( {(1,0):F(11,23)}, ),( {(0,1):F(13,29)}, ))
    modes=((1,2,True,True),(2,1,False,True),(1,2,True,False),(2,2,False,False))
    for full,single,cf,cs in modes:
        def reduce(table):
            out=[]
            for row in table:
                result=[]
                for item in row:
                    v=defaultdict(F)
                    for (f,s),x in item.items():
                        if cf: f%=full+1
                        if cs: s%=single+1
                        if f<=full and s<=single: v[f,s]+=x
                    result.append(dict(v))
                out.append(tuple(result))
            return tuple(out)
        def exact_product(left,right):
            out=[]
            for row in left:
                result=[]
                for j in range(len(right[0])):
                    item=defaultdict(F)
                    for k,p in enumerate(row):
                        for (f,s),x in p.items():
                            for (ff,ss),y in right[k][j].items(): item[f+ff,s+ss]+=x*y
                    result.append(dict(item))
                out.append(tuple(result))
            return reduce(out)
        arith=IntegerArithmetic(full,single,cf,cs,bits=8,guard=guard)
        rr,tt=arith.pack_matrix(a),arith.pack_matrix(b)
        er,et=reduce(a),reduce(b)
        for _ in range(6):
            guard(); tt=arith.multiply(rr,tt); rr=arith.multiply(rr,rr)
            et=exact_product(er,et); er=exact_product(er,er)
            for actual,expected in ((rr,er),(tt,et)):
                for s in range(single+1):
                    for f in range(full+1):
                        values=arith.extract(actual,f,s)
                        for i,row in enumerate(values):
                            for j,value in enumerate(row):
                                assert value>=expected[i][j].get((f,s),F(0)); checks+=1
                mass=arith.mass(actual)
                for i,row in enumerate(mass):
                    for j,value in enumerate(row):
                        assert value>=sum(expected[i][j].values(),F(0)); checks+=1
        assert arith.products==12 and arith.rounded_coefficients>0
        # Overflow both axes: fold exactly where requested, otherwise discard.
        left=arith.pack_matrix((( {(full,single):F(1)}, ),))
        right=arith.pack_matrix((( {(1,1):F(1)}, ),))
        carried=arith.multiply(left,right)
        for s in range(single+1):
            for f in range(full+1):
                wanted=F(int(cf and cs and f==s==0))
                assert arith.extract(carried,f,s)[0][0]==wanted; checks+=1
        zero=arith.pack_matrix((( {}, ),))
        zero=arith.multiply(zero,zero)
        assert zero.exponent==0 and arith.mass(zero)==((F(0),),); checks+=1
    # Constant noncommuting matrices directly test terminal order at64 factors.
    arith=IntegerArithmetic(0,0,bits=128,guard=guard)
    caps={(0,0):tuple(map(F,(1,2,3,4,5,7)))}
    rr,tt=arith.pack(caps)
    for _ in range(6): tt=arith.multiply(rr,tt); rr=arith.multiply(rr,rr)
    exact_r=((F(1),F(2)),(F(3),F(4))); exact_t=(F(5),F(7))
    for _ in range(63):
        exact_t=(exact_r[0][0]*exact_t[0]+exact_r[0][1]*exact_t[1],
                 exact_r[1][0]*exact_t[0]+exact_r[1][1]*exact_t[1])
    for i in range(2): assert arith.extract(tt,0,0)[i][0]>=exact_t[i]; checks+=1
    return {"exact_fraction_comparisons":checks,"six_doubling_layouts":len(modes),
        "coarse_mantissa_bits_for_rounding_tests":8,"terminal_order":"H^63 J",
        "carry_zero_and_noncommuting_controls":True}


def peak_memory_reader():
    if sys.platform=="win32":
        import ctypes
        from ctypes import wintypes
        class Counters(ctypes.Structure):
            _fields_=[("cb",wintypes.DWORD),("PageFaultCount",wintypes.DWORD),
                *[(x,ctypes.c_size_t) for x in ("PeakWorkingSetSize","WorkingSetSize",
                "QuotaPeakPagedPoolUsage","QuotaPagedPoolUsage","QuotaPeakNonPagedPoolUsage",
                "QuotaNonPagedPoolUsage","PagefileUsage","PeakPagefileUsage")]]
        process=ctypes.windll.kernel32.GetCurrentProcess; process.restype=wintypes.HANDLE
        info=ctypes.windll.psapi.GetProcessMemoryInfo
        info.argtypes=[wintypes.HANDLE,ctypes.POINTER(Counters),wintypes.DWORD]
        def peak():
            counters=Counters(); counters.cb=ctypes.sizeof(counters)
            if not info(process(),ctypes.byref(counters),counters.cb): raise OSError("memory guard unavailable")
            return counters.PeakWorkingSetSize
    else:
        import resource
        def peak():
            n=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            return n if sys.platform=="darwin" else 1024*n
    return peak


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); peak=peak_memory_reader()
    def guard():
        if time.monotonic()-started>args.seconds: raise TimeoutError("integer regional test deadline; no retry")
        if peak()>MEMORY: raise MemoryError("512MiB integer regional test budget exceeded")
    tests=exact_tests(guard); guard()
    print(json.dumps({"scope":"exact small arithmetic tests only; no regional certificate",
        "tests":tests,"seconds":time.monotonic()-started,"peak_process_bytes":peak(),
        "generated_bytes":0,"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
