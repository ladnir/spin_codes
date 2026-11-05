"""Whole-word Minkowski for a labelled singleton RS row.

Each original column is either zero or nonzero. Keep the >=9 active
columns restriction BEFORE taking a whole-word Lr norm. This is a
complementary positive bound, not necessarily smaller than entrywise
block Minkowski: summing state paths before taking a root can cost slack.
No regional factor, RS normalization, source count, or row-set sum here.
"""
import numpy as np


def whole_word_roots(blocks, groups=18, *, minimum_activity=9, q=256,
                     check_time=None):
    """Per-A log2 upper norm, summed over allowed original activity words.

    Eight blocks have shape (number_of_A,4,2,2). Branch order is
    (00,01,10,11), where the left column is the lower-numbered column.
    J_block entries contain the r-th moment of that fixed branch.
    Both q and groups are positive integer parameters of the comparison.
    """
    if (not 1 <= len(blocks) <= 8 or type(groups) is not int or groups < 1
            or type(q) is not int or q < 2 or type(minimum_activity) is not int
            or not 0 <= minimum_activity <= 2*len(blocks)):
        raise ValueError('Bounded blocks, valid field size, exponent and activity required')
    matrices = [np.asarray(block, dtype=float) for block in blocks]
    shape = matrices[0].shape
    if (len(shape) != 4 or shape[1:] != (4,2,2)
            or any(block.shape != shape or np.any(np.isnan(block))
                   or np.any(np.isposinf(block)) for block in matrices)):
        raise ValueError('Compatible labelled-A x four-branch x two-state blocks required')
    masks = np.asarray([mask for mask in range(1 << (2*len(blocks)))
                        if mask.bit_count() >= minimum_activity], dtype=np.uint32)
    weights = np.asarray([int(mask).bit_count() for mask in masks]) * np.log2(q-1.)
    indices = [2*((masks >> (2*b)) & 1) + ((masks >> (2*b+1)) & 1)
               for b in range(len(blocks))]
    answer = np.empty(shape[0])
    for source in range(shape[0]):
        state = np.full((len(masks),2), -np.inf)
        state[:,0] = 0.
        for block, index in zip(matrices, indices):
            if check_time is not None:
                check_time()
            current = block[source,index]
            state = np.logaddexp2(state[:,0,None]+current[:,0,:],
                                 state[:,1,None]+current[:,1,:])
        norms = np.logaddexp2.reduce(state,axis=1)/groups
        answer[source] = np.logaddexp2.reduce(weights+norms)
    return answer
