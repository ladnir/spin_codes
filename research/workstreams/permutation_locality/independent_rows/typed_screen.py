"""Binary64 proposal-only placement of two persistent group types.

This is a numerical screening analogue of typed_placement. It returns
diagnostic matrices, never outward upper bounds or proof certificates.
Every promising witness must be replayed by the exact/outward driver.

The normalized hypergeometric recurrence matches typed_placement, including
selected groups whose reference packets happen to be zero. For each local
type count, NumPy batches all prefix-count matrix products. Matrices need
not commute: each previous-prefix operator multiplies its epoch operator
on the right. Two count-grid buffers and one matrix-product scratch buffer
are reused across epochs.
"""
import numpy as np


RESULT_SCOPE='Binary64 proposal only; no outward certificate.'


def _as_array(value):
    if hasattr(value,'nrows') and hasattr(value,'ncols'):
        value=[[float(value[i,j]) for j in range(value.ncols())]
               for i in range(value.nrows())]
    return np.asarray(value,dtype=np.float64)


def _inputs(operators, counts, epochs, slots):
    if (not isinstance(epochs,int) or epochs<0
            or not isinstance(slots,int) or slots<1
            or len(counts)!=2 or any(not isinstance(x,int) or x<0 for x in counts)
            or sum(counts)>epochs*slots or not operators):
        raise ValueError('invalid epochs, slots, type counts, or empty operators')
    result={}
    size=None
    for key,value in operators.items():
        if (len(key)!=2 or any(not isinstance(x,int) or x<0 for x in key)
                or sum(key)>slots):
            raise ValueError('invalid local type counts')
        value=_as_array(value)
        if (value.ndim!=2 or value.shape[0]<1 or value.shape[0]!=value.shape[1]
                or not np.isfinite(value).all() or (value<0).any()):
            raise ValueError('finite nonnegative square epoch matrices required')
        if size is None:
            size=value.shape[0]
        elif value.shape!=(size,size):
            raise ValueError('epoch matrices must have matching dimensions')
        result[key]=value
    return result,size


def _log_choose_table(positions, slots):
    """Stable small-degree log binomials without gamma subtraction."""
    table=np.full((slots+1,positions+1),-np.inf,dtype=np.float64)
    table[0,:]=0.
    population=np.arange(positions+1,dtype=np.float64)
    for chosen in range(1,slots+1):
        table[chosen,chosen:]=(table[chosen-1,chosen:]
                               +np.log((population[chosen:]-chosen+1)/chosen))
    return table


def screen_placement(operators, counts, *, epochs=64, slots=32, return_grid=False):
    """Return a proposal-only NumPy matrix, or the full bounded count grid.

    operators[(j1,j2)] may be a NumPy, Arb, or exact rational matrix; it is
    converted to binary64. A returned scalar entry is not a certified bound.
    return_grid=True returns shape (n1+1,n2+1,d,d). The default returns the
    matrix for counts=(n1,n2) and prunes prefixes unable to reach that pair.
    """
    counts=tuple(counts)
    operators,size=_inputs(operators,counts,epochs,slots)
    n1,n2=counts
    shape=(n1+1,n2+1,size,size)
    current=np.zeros(shape,dtype=np.float64)
    current[0,0]=np.eye(size)
    if not epochs:
        return current if return_grid else current[0,0].copy()
    updated=np.empty_like(current)
    scratch=np.empty_like(current)
    totals=np.arange(n1+1)[:,None]+np.arange(n2+1)[None,:]
    log_choose=_log_choose_table(epochs*slots,slots)
    terms=[]
    for j1 in range(min(n1,slots)+1):
        for j2 in range(min(n2,slots-j1)+1):
            if not return_grid and n1+n2-j1-j2>(epochs-1)*slots:
                continue
            if (j1,j2) not in operators:
                raise ValueError(f'missing epoch operator for type counts {(j1,j2)}')
            terms.append((j1,j2,operators[j1,j2],slots-j1-j2))
    for e in range(1,epochs+1):
        positions=e*slots
        valid=totals<=positions
        if not return_grid:
            valid&=n1+n2-totals<=(epochs-e)*slots
        normalization=np.zeros((n1+1,n2+1),dtype=np.float64)
        updated.fill(0.)
        log_denominator=log_choose[slots,positions]
        for j1,j2,operator,empty_slots in terms:
            target_valid=valid[j1:,j2:]
            if not target_valid.any():
                continue
            empty=np.maximum(positions-totals[j1:,j2:],0)
            log_weight=(log_choose[j1,j1:n1+1,None]
                        +log_choose[j2,None,j2:n2+1]
                        +log_choose[empty_slots,empty]-log_denominator)
            log_weight=np.where(target_valid,log_weight,-np.inf)
            weights=np.exp(log_weight)
            normalization[j1:,j2:]+=weights
            a,b=weights.shape
            block=scratch[:a,:b]
            np.matmul(current[:a,:b],operator,out=block)
            np.multiply(block,weights[:,:,None,None],out=block)
            np.add(updated[j1:,j2:],block,out=updated[j1:,j2:])
        # A diagnostic guard against indexing/normalization mistakes, not
        # an interval-arithmetic certificate for floating computations.
        if not np.allclose(normalization[valid],1.,rtol=0.,atol=1e-10):
            raise ArithmeticError('proposal hypergeometric normalization failed')
        if not np.isfinite(updated).all():
            raise ArithmeticError('proposal placement overflowed or produced a nonfinite matrix')
        current,updated=updated,current
    return current if return_grid else current[n1,n2].copy()
