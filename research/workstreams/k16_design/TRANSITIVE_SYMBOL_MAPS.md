# Cheaper Symbol Maps with the Same First-Moment Bound

The larger-outer proof uses each randomized nonzero symbol's distribution,
not the full distribution of its binary matrix. This permits two exact
changes to the outer. Neither change requires a heuristic distance estimate.
The remaining routing and inner distributions stay unchanged.

## Eight Nibble Rows Instead of Four Byte Rows

Use eight parallel RS[16,8] rows over GF(16), with common evaluation points.
An aligned symbol contains eight nibbles, or 32 binary coordinates. Each
group therefore maps 256 input bits to 512 output bits and supplies 128
four-bit packets. These are the same group dimensions as four parallel
RS[16,8] rows over GF(256).

Identify GF(16)^8 with a degree-eight extension field E through a fixed
basis. The common GF(16) generator acts coordinatewise, so the aligned
code is an E-linear MDS[16,8] code with |E| = 2^32. For each fixed symbol
support S, shortening has E-dimension max(0, |S|-8). Inclusion-exclusion
therefore determines its exact-support count solely from |S| and |E|.
Four parallel GF(256) rows give the same counts through a degree-four
extension. No specific field isomorphism needs to be evaluated by the encoder.

Apply independent uniform nonzero-symbol randomizers at the 16 positions.
For any fixed support S, every output with precisely that support has the
same expected multiplicity. The exact-support counts above show that the
two outers have identical expected measures on all binary output words.
Consequently they have the same pointwise uniform majorant and the same
expected packet-support enumerator. The existing larger-outer first-moment
bounds can use either implementation after their construction metadata is
updated; this does not transfer a certificate for a different inner state.

The engineering attraction is concrete: the eight-row implementation can
reuse the existing factored nibble RS circuit instead of developing a new
byte RS circuit. This statement predicts structure reuse, not a timing.

## Field Multiplication Instead of a Uniform GL32 Matrix

Fix an arbitrary binary identification of 32-bit strings with a field F
of size 2^32. At each symbol position, independently sample c from F minus
zero and apply the binary linear map x -> c*x. Sampling occurs once during
setup, before the realized code is used for messages.

For every fixed x != 0, multiplication by x is a bijection on F minus zero.
Thus c*x is uniform on all nonzero 32-bit strings. For x = 0, the output
is zero. These are exactly the two fixed-input laws of a uniform GL32
matrix. Independence of the scalars supplies the same joint law across
symbol positions for each fixed message. Summing over messages preserves
the expected outer measure and the first-moment distance bound.

The scalar family is not uniform in GL32. In particular, this argument
does not preserve correlations between different messages or automatically
transfer a second-moment proof. Those correlations are not used here.
The scalar multiplication basis need not agree with the abstract extension
basis used to count the RS words: both constructions randomize each nonzero
binary symbol over the same set.

## A Nine-Map Binary-Adjoint Circuit

The prototype represents F using a quadratic tower over the byte field
GF(2)[z]/(z^8+z^4+z^3+z+1). Let mu = 0x20 and define

    u^2 + u + mu = 0,
    v^2 + v + mu*u = 0.

The four coordinate bytes have basis (1,u,v,uv), least significant first.
The absolute trace of mu in GF(256) is one. Since the relative trace of u
over GF(256) is one, the absolute trace of mu*u in GF(65536) is also one.
Both quadratic polynomials are therefore irreducible.

For a quadratic extension w^2+w+eta=0, multiplying x0+x1*w by a fixed
c0+c1*w can be evaluated as

    p0 = c0*x0,
    p1 = (eta*c1)*x1,
    p2 = (c0+c1)*(x0+x1),
    output = (p0+p1, p0+p2).

All operations are in characteristic two. Reduction is folded into the
setup constant eta*c1. Recursing through the two quadratic extensions
gives nine constant byte maps and fifteen byte-vector XORs. No separate
reduction map is needed in the encoding path.

Transposed encoding uses the binary adjoint of this circuit, not ordinary
multiplication by the scalar. Each constant byte map is replaced by its
8-by-8 binary transpose, and circuit edges are reversed. The implementation
uses nine GFNI affine instructions and fifteen vector XORs for each payload
half, with 72 bytes of precomputed coefficients per symbol. A dense GL32
uses 128 coefficient bytes and sixteen affine instructions per half.
These operation counts do not establish an end-to-end speedup.

`field_symbol_randomizer.py` provides independent scalar arithmetic and
adjoint checks. `implementation/rs16x8/Tower32Randomizer.h` supplies setup
helpers and the explicit SIMD circuit. The caller must independently
sample a uniform nonzero 32-bit scalar at each symbol; rejecting zero from
uniform 32-bit words has no modulo bias.

## A Larger Group Without Longer Nibble RS Rows

Sixteen parallel GF(16) RS[16,8] rows give a 512-to-1,024-bit group.
Apply an independent uniform nonzero GF(2^64) multiplier to each aligned
sixteen-nibble symbol. Then split into 256 four-bit packets, shuffle the
packets independently within each group, and shuffle every region independently.
At K = 2^20 there are 2,048 groups, the same geometry as the retained
four-row GF(256) RS[32,16] proof anchor. The inner remains unchanged.

The aligned outer is an MDS[16,8] code over an effective alphabet of size
2^64. On any exact support S of size h >= 9, shortening leaves dimension
h-8. An information set of that dimension injects fully supported words
into h-8 nonzero symbols, so their number is at most (2^64-1)^(h-8).
Independent symbol randomizers consequently give the pointwise majorant

    beta64 = 2^1024 / (2^64-1)^8.

The retained GF(256) RS[32,16] construction uses

    beta32 = 2^1024 / (2^32-1)^16,
    beta64 / beta32 = ((2^32-1)/(2^32+1))^8 < 1.

Thus every uniform-majorant occupancy bound proved with beta32 remains
valid for this new outer, at the same length, inner, and routing geometry.
The retained larger-outer certificate uses this majorant for q1 and q2
as well as its transferred tail, so this comparison covers all occupancies.
It does not transfer a bound based on a different exact outer-shell
enumerator or change the certified inner state size.

This alternative does not require RS[32,16] over GF(16): the rows remain
length 16, and the group grows by adding parallel rows.

The field implementation adds a third quadratic extension to the preceding
tower, with w^2+w+theta=0 and theta=mu*u*v=0x20000000. The absolute trace
of theta is one by the same relative-trace argument. The eight byte
coordinates are (1,u,v,uv,w,uw,vw,uvw). Recursing the constant-multiplier
circuit gives 27 byte GFNI maps and 57 vector XORs per payload half,
with 216 coefficient bytes per symbol. `Tower64Randomizer.h` retains this
explicit binary-adjoint circuit. Runtime and full-encoder validation are
separate from the algebraic checks.

## Moving the Adjoint into the Sampled Family

An alternative samples the forward symbol map as M_c^T, where M_c is
ordinary multiplication by an independently uniform nonzero field element c.
The transposed encoder then evaluates M_c directly. This family also sends
each fixed nonzero binary vector uniformly to all nonzero vectors.

To see this, fix x != 0. If M_c^T*x = M_d^T*x for distinct c and d, then
M_(c+d)^T*x = 0. Multiplication by c+d is invertible, so its transpose
has trivial kernel, a contradiction. Thus c -> M_c^T*x is injective and
hence bijective on the nonzero field elements. Independent sampling at each
symbol again gives the fixed-message law used by the first-moment proof.

`Tower32ByteRandomizer.h` implements the resulting transposed symbol map
with nine GFNI byte multiplications and fifteen vector XORs per payload
half. Its coefficient table stores nine bytes rather than nine 64-bit
affine matrices. It has the same asymptotic arithmetic count but different
instruction and memory costs, which require an end-to-end comparison.
Its scalar oracle stores the rows of ordinary multiplication, not the
rows of the adjoint multiplier used by `Tower32Randomizer.h`.
