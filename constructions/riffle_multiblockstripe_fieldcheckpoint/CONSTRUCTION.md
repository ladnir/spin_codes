# Construction

Fix (B=256) outer coordinates, (M=8192) outer blocks, and a divisor
(q\mid32).  Write (S:=M/q).  The construction groups the outer blocks into
(S) groups of (q) blocks.

Each outer block first receives an independent permutation of its (B)
coordinates.  For coordinate (c\in[B]), block lane (j\in[q]), and group
(u\in[S]), the route assigns the bit to lane slice

\[
 r:=j+c+\delta_{u,c}\pmod q.
\]

The route then applies an independent permutation of the (S) group positions
inside each pair ((c,r)).  The output order is coordinate-major and then
lane-major.  Each outer block therefore contributes one candidate position to
every 8192-bit transposed region.

Two shift schedules define distinct variants.

- **Riffle MultiBlockStripe-(q) FieldCheckpoint** samples one
  \(\delta_u\in[q]\) per group and sets \(\delta_{u,c}:=\delta_u\).
- **Riffle MultiBlockStripeFresh-(q) FieldCheckpoint** samples every
  \(\delta_{u,c}\in[q]\) independently.

Both variants use the same precomputed routing interface.  Fresh shifts change
setup and schedule generation, but they add no operation to the encoder.  At
(q=32), the fresh shifts consume (256\cdot256\cdot5=327{,}680) setup bits,
or 40 KiB.

The remaining stages match Riffle ExactPerm FieldCheckpoint v1.  The outer code
is a binary ([256,128,\mathord\ge38]) code.  FieldCheckpoint has 64 state bits
and a 256-bit epoch.

