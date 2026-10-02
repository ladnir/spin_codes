"""Explicit t64/s10 expansion from three GF8 bilinear product coordinates.

Evaluation points z are the integers0..63 in low-bit-first coordinates.
X=(z0,z2,z4), Y=(z1,z3,z5), and GF8 uses modulus x^3+x+1.
Rows are1,z0,..,z5, followed by the three low-bit-first coordinates of X*Y.
Feedback is exactly the transpose. The inner samples fresh uniform GL10
at every physical64-bit step, starts at zero, and has no reset or flush.
"""
from collections import Counter
import hashlib
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import screen_t128 as screen


def multiply8(a,b):
    result=0
    for _ in range(3):
        if b&1:result^=a
        b>>=1
        a=(a<<1)^(0xb if a&4 else 0)
    return result


def product_at(z):
    x=sum(((z>>(2*i))&1)<<i for i in range(3))
    y=sum(((z>>(2*i+1))&1)<<i for i in range(3))
    return multiply8(x,y)


def prepare():
    algebra=screen.q1.kernel_t64.s16_maps
    rows=[(1<<64)-1]+[sum(((z>>i)&1)<<z for z in range(64)) for i in range(6)]
    rows += [sum(((product_at(z)>>i)&1)<<z for z in range(64)) for i in range(3)]
    rows=tuple(rows)
    columns=tuple(sum(((row>>z)&1)<<i for i,row in enumerate(rows)) for z in range(64))
    polar_ranks=[]
    for selector in range(1,8):
        q=lambda z:(product_at(z)&selector).bit_count()%2
        polar=[sum((q((1<<i)^(1<<j))^q(1<<i)^q(1<<j)^q(0))<<j for j in range(6)) for i in range(6)]
        polar_ranks.append(algebra.binary_rank(polar))
    packet_ranks=[algebra.binary_rank(columns[i:i+4]) for i in range(0,64,4)]
    if (algebra.binary_rank(rows)!=10 or algebra.binary_rank(columns)!=10 or
            polar_ranks!=[6]*7 or packet_ranks!=[4]*16 or
            any((a&b).bit_count()&1 for a in rows for b in rows)):
        raise ArithmeticError('rank, nondegenerate polar form, packet rank, or CA check failed')
    images=tuple(algebra.images_from_rows(rows))
    spectrum=dict(sorted(Counter(x.bit_count() for x in images).items()))
    expected={0:1,28:448,32:126,36:448,64:1}
    if spectrum!=expected or len(set(images))!=1024:
        raise ArithmeticError('complete injective state census or expected spectrum failed')
    physical=screen.q1.kernel_t64.kernel_maps.prepare_maps(images,columns,bits=10,
        distribution='uniform_gl',birth_density='capped')
    record=dict(schema='bilinear-gf8-t64-s10-fixed-map-1',physical_t=64,state_bits=10,
        field_modulus='x^3+x+1',X_coordinates=[0,2,4],Y_coordinates=[1,3,5],
        expansion_rows_hex=list(map(hex,rows)),feedback_columns=list(columns),
        expansion_rank=10,feedback_rank=10,feedback_definition='C=A^T',
        feedback_times_expansion_zero=True,packet_ranks=packet_ranks,
        nonzero_quadratic_polar_ranks=polar_ranks,state_count=1024,
        full_state_census=True,expansion_spectrum={str(k):v for k,v in spectrum.items()},
        minimum_expansion_weight=28,map_sha256=physical['map_sha256'],
        state_update='independent uniform GL10 per physical64-bit step',
        macro=dict(physical_steps=2,physical_step_bits=64,macro_step_bits=128,
                   state_continuity='retained_between_halves'))
    return screen.q1.kernel_t64.wrap(physical),record


def source_pins():
    pins=screen.q1.source_snapshot()
    for path in (Path(screen.__file__).resolve(),Path(__file__).resolve(),HERE/'screen.py'):
        pins[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    return pins
