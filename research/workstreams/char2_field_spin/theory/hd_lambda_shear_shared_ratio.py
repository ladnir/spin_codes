"""A bounded FLOAT common-atom majorant for the scheduled column shear.

Fix two projective pre-shear vectors U,V.  Their relative scale rho is
uniform in F_q^*: uniform column inputs, or the full independent D1 law,
permit this conditioning.  Only ONE common pair gauge may be removed.
Even physical rows use (U_i+rho V_i,rho V_i); odd rows use
(U_i,U_i+rho V_i).  The final independent D2 remains after the shear.

At output marks in [0,1], let M_OR retain every overlapping target and
M_XOR cancel every overlapping target.  The exact rho-average obeys

    E_rho M <= (q-2)/(q-1) M_OR + 1/(q-1) M_XOR.

To prove this, partition overlapping rows by their cancellation ratio.
Merge two disjoint classes C,D into one exceptional class.  Replacing
their two terms by a union-cancel term and an OR term increases the sum
by unchanged_product*(1-product(t_C))*(1-product(t_D)) >= 0.
Iterate.  This proof includes zero marks and needs no division by them.
It is NOT coefficientwise support domination or a density/Renyi bound.

Use the same positive matrix atom measure at all positions, including
its artificial empty transitions.  For fixed block boundary states let
J_OR,J_XOR be the r-th support-monomial moments.  Minkowski gives the
entrywise norm upper (q-2)/(q-1)*J_OR^(1/r)+1/(q-1)*J_XOR^(1/r).
The two coefficients stay OUTSIDE the roots; atom masses are unpowered.
Dropping binomial cross terms inside an r-th moment would be invalid.

For fixed original A, single-column histograms n_j(w,S) are exact.
They retain original input weight and physical Q_j output labels.  The
product counts on consecutive columns may use the ratio inequality:
scaling either nonzero column preserves its input weight and support.
No extra q-1 factor is charged.  A's label set stays fixed throughout
the sixteen columns; relation kernels are not held fixed across them.

pair_coefficients returns coefficients of the resulting branch-root
polynomial.  The caller raises EACH entry to r, multiplies eight matrix
polynomials chronologically, and retains degree >=9*|A|*r.  The RS
factor 255^(-8|A|), row-set root sum, and common prefix are NOT included.
There is no gate, regional provider, native operation, or certificate.
"""
from dataclasses import dataclass
from hashlib import sha256
from math import log2
from pathlib import Path

import numpy as np

import hd_lambda_scheduled_relations as scheduled
from hd_pair_identity_kernel import mask_table_from_entries, branch_log_matrices
from hd_pair_state_kernel import common_moment_atoms


DIRECTORY = Path(__file__).resolve().parent
SCHEDULE = 'input-inverse-output-xor'
MAX_HISTOGRAM = 529
MAX_CHUNK = 8192
PINS = {
    'hd_lambda_scheduled_relations.py': 'eeccb542a259b60f6e37dc6ec8aa8b0b724da78483a1bf0e0b5f084727739e0e',
    'hd_lambda_outer.py': 'f0c063843e295a160a1186c6c084a07dd1c2a3745005ed68057fde59913c99e7',
    'hd_pair_scheduled_outer.py': '64dc50bf18171d371ae16f4484087e9d42d89cc7c9b55a7576be55b1251e4741',
    'hd_quaternary_outer_relations.py': 'd8be4e7a1ab77b98c76695f13c1bc9a9cf867df2aae6408ac3e08c6fc12346c0',
    'hd_outer_pair_layers.py': '017c5ba20747ae8e25acf02393e9702a03b086038e1814ac2424cdfe49d66f60',
    'hd_binary_syndrome_outer.py': 'b1dc2dcb6e212a28526f7c7dc7f7a55999be7e6f5f53370f78768d705f2bda96',
    'hd_pair_identity_kernel.py': '93f3a9053a72bed39db6b00a32d2f2b7b060a35ce19ce814dd856d096a796961',
    'hd_pair_state_kernel.py': 'b28142e19e51a906ae742128325d1e7557dfb425624c6da45778a2a75753e195',
    'hd_post_refresh_quiet_rows.py': '75e38819d225e9771ade1ed90509b7d4582a2637b24b6522f5b62b0200eb49e7',
}


def dependencies():
    result = dict(PINS)
    for name, digest in result.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen dependency changed: '+name)
    result[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


def _logs(values):
    values = np.asarray(values, dtype=float)
    if np.any(np.isnan(values)) or np.any(np.isposinf(values)):
        raise ValueError('Finite log2 values or exact zero(-inf) required')
    return values


@dataclass(frozen=True)
class RatioKernel:
    """A prepared common measure; create with prepare_atoms/prepare_measure."""
    table: np.ndarray
    groups: int
    rows: int


def prepare_atoms(log_atoms, log_matrices, groups, *, rows=16):
    """Prepare r-th mask moments once, rejecting any scalar atom above one."""
    atoms, matrices = _logs(log_atoms), _logs(log_matrices)
    if (atoms.ndim != 1 or not 1 <= len(atoms) <= 32
            or matrices.shape != (len(atoms), 2, 2) or np.any(atoms > 0)
            or type(groups) is not int or not 1 <= groups <= 18
            or type(rows) is not int or not 1 <= rows <= 16):
        raise ValueError('Atoms in[0,1], positive matrix masses, groups1..18, rows1..16 required')
    zero = np.logaddexp2.reduce(matrices, axis=0)
    nonzero = np.logaddexp2.reduce(matrices+groups*atoms[:, None, None], axis=0)
    table = mask_table_from_entries(zero, nonzero, rows=rows)
    table.setflags(write=False)
    return RatioKernel(table, groups, rows)


def prepare_measure(measure, groups=18):
    """Use the unchanged five-atom matrix measure; do not normalize its mass."""
    return prepare_atoms(*common_moment_atoms(measure), groups, rows=16)


def output_masks(first, second, *, rows=16):
    """Return (OR-left,OR-right,XOR-left,XOR-right) physical support masks."""
    first, second = np.asarray(first), np.asarray(second)
    if (type(rows) is not int or not 1 <= rows <= 16
            or first.shape != second.shape
            or not np.issubdtype(first.dtype, np.integer)
            or not np.issubdtype(second.dtype, np.integer)
            or np.any(first < 0) or np.any(second < 0)
            or np.any(first >= 1 << rows) or np.any(second >= 1 << rows)):
        raise ValueError('Equal-shaped integer support masks in the declared rows required')
    even = sum(1 << i for i in range(0, rows, 2))
    odd = ((1 << rows)-1) ^ even
    union, difference = first | second, first ^ second
    return ((union & even) | (first & odd),
            (second & even) | (union & odd),
            (difference & even) | (first & odd),
            (second & even) | (difference & odd))


def mixture_log_norms(first, second, prepared, *, q=256):
    """Batched log2 norm upper, with ratio probabilities outside r-th roots.

    This is not log J.  Do not divide its result by r a second time.
    Matrix entries already sum all atom/state paths within both columns.
    """
    if (not isinstance(prepared, RatioKernel) or type(q) is not int
            or q not in (4, 16, 256)):
        raise ValueError('A prepared monotone kernel and GF4-extension field required')
    masks = output_masks(first, second, rows=prepared.rows)
    or_moment = branch_log_matrices(masks[0], masks[1], prepared.table)
    xor_moment = branch_log_matrices(masks[2], masks[3], prepared.table)
    return np.logaddexp2(log2(q-2)-log2(q-1)+or_moment/prepared.groups,
                        -log2(q-1)+xor_moment/prepared.groups)


def column_histogram(source, column, *, q=256, schedule=SCHEDULE):
    """Exact (original input weight, physical output mask, Python-int count)."""
    counts = scheduled.joint_counts(source, schedule, column, q=q)
    histogram = tuple((weight, mask, count) for (weight, mask), count in counts.items())
    if sum(row[2] for row in histogram) != q**len(source):
        raise ArithmeticError('Single-column raw mass changed')
    return histogram


def _histogram(histogram, rows):
    histogram = tuple(histogram)
    if (not 1 <= len(histogram) <= MAX_HISTOGRAM
            or any(len(row) != 3 or any(type(x) is not int for x in row)
                   or not 0 <= row[0] <= 4 or not 0 <= row[1] < 1 << rows
                   or row[2] <= 0 for row in histogram)
            or len({(row[0], row[1]) for row in histogram}) != len(histogram)):
        raise ValueError('At most529 distinct positive exact column histogram rows required')
    return histogram


def pair_coefficients(first_histogram, second_histogram, prepared, *, q=256,
                      chunk_size=8192, check_time=None):
    """Log polynomial coefficients, shape(2a+1,2,2), ready for entrywise r-power.

    Exact integer counts are retained until their individual log conversion;
    products are formed by adding logs, not overflowing a uint64 array.
    At most529**2 support pairs are processed in bounded vectorized chunks.
    """
    if (not isinstance(prepared, RatioKernel) or type(chunk_size) is not int
            or not 1 <= chunk_size <= MAX_CHUNK):
        raise ValueError('Prepared kernel and chunk size1..8192 required')
    first = _histogram(first_histogram, prepared.rows)
    second = _histogram(second_histogram, prepared.rows)
    wa = np.array([row[0] for row in first], dtype=np.int8)
    wb = np.array([row[0] for row in second], dtype=np.int8)
    sa = np.array([row[1] for row in first], dtype=np.uint16)
    sb = np.array([row[1] for row in second], dtype=np.uint16)
    na = np.array([log2(row[2]) for row in first])
    nb = np.array([log2(row[2]) for row in second])
    maximum = int(wa.max()+wb.max())
    result = np.full((maximum+1, 2, 2), -np.inf)
    width = len(second)
    for start in range(0, len(first)*width, chunk_size):
        if check_time is not None:
            check_time()
        indices = np.arange(start, min(start+chunk_size, len(first)*width))
        left, right = indices//width, indices%width
        terms = mixture_log_norms(sa[left], sb[right], prepared, q=q)
        terms += (na[left]+nb[right])[:, None, None]
        weights = wa[left]+wb[right]
        for degree in range(maximum+1):
            selected = terms[weights == degree]
            if len(selected):
                result[degree] = np.logaddexp2(result[degree],
                                             np.logaddexp2.reduce(selected, axis=0))
    return result


def pair_polynomial(source, column, prepared, *, q=256, schedule=SCHEDULE,
                    chunk_size=8192, check_time=None):
    """One original labelled A and consecutive physical columns j,j+1."""
    if (type(column) is not int or column not in range(0, 16, 2)
            or not isinstance(prepared, RatioKernel) or prepared.rows != 16):
        raise ValueError('An even column0..14 and a sixteen-row prepared kernel required')
    return pair_coefficients(column_histogram(source, column, q=q, schedule=schedule),
        column_histogram(source, column+1, q=q, schedule=schedule), prepared,
        q=q, chunk_size=chunk_size, check_time=check_time)


def metadata():
    return dict(arithmetic='FLOAT-log2-not-directed', supported_arities=[1,2,3,4],
        ratio_shared_across_all_rows=True, atom_domain='[0,1]',
        ratio_coefficients_outside_roots=True, atom_masses_powered=False,
        pointwise_common_atom_majorant=True, coefficientwise_support_domination=False,
        density_or_renyi_upper=False, original_A_fixed_across_columns=True,
        physical_even_target='left', physical_odd_target='right',
        full_independent_D1_comparison=True, permitted_pair_gauges=1,
        final_D2_after_shear=True, rs_normalization_included=False,
        source_count_included=False, regional_provider_runs=0,
        production_gate_runs=0, whole_code_certificate=False)
