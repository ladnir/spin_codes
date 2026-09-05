# Goal 10: close the all-one endpoint family

## Question

Goal 09 reduced the finite-support endpoint contribution from

\[
2^{159.773}
\quad\text{to}\quad
2^{40.744}.
\]

That calculation still treated every packet in the all-one BCH block as a
free packet value. Can the exact packet structure close this subcase?

The output threshold is

\[
D=188{,}766=\lfloor 0.09L\rfloor,
\qquad L=2{,}097{,}408.
\]

## Fixed-packet operator

Every all-one BCH block contains 32 packets of value `1111`. Local bit
permutation leaves these packets unchanged.

First consider a finite support from Goal 09. One BCH block is all ones. The
other two blocks are complementary and have total binary weight 128. After
the global packet permutation, the active superblock has this description:

- 32 fixed packets have value `1111`;
- 64 variable packets contain 128 one-bits in total;
- the fixed and variable packet labels have a uniform random interleaving.

The new dynamic program retains

\[
(f,H,q),
\]

where \(f\) counts placed fixed packets, \(H\) is the active packet support,
and \(q\in\mathbb F_2^4\) is the accumulator state. The program ends at
\(f=32\). It divides the path sum by \(\binom{96}{32}\), the number of
fixed-packet interleavings.

A coefficient tilt enforces total weight 128 only on the 64 variable
packets. Thus, for every \(x>0\), the operator uses

\[
[x^{128}]F(x)\le x^{-128}F(x).
\tag{1}
\]

Goal 09 instead applied the coefficient tilt to all 96 packets at total
weight 256. The new operator removes that false family while retaining the
same 16 accumulator states and the same global gap bound.

The implementation passes two independent checks. Direct enumeration agrees
on an instance with two fixed and three variable packets. With neutral
weights, the full operator has total mass

\[
\binom{96}{32}16^{64}.
\]

## Finite-support result

Transferring the complete Goal 09 histogram through the fixed-packet
operator gives

\[
\mathbb E Z_D^{\mathrm{finite,end}}
\le 2^{-44.5527}.
\tag{2}
\]

Preserving the all-one block therefore gains 85.297 bits over the Goal 09
transfer. Equation (2) closes every endpoint word on a support that excludes
the second parity position.

## Supports containing the second parity position

Let a support contain the second parity position and two finite projective
points with coefficients \(a\ne b\). The two finite field values are equal.
There are two endpoint normalizations.

1. A finite coordinate is all ones. Both finite BCH blocks are then all
   ones, and the parity value is \(e(a+b)\).
2. The second parity coordinate is all ones. Both finite values are
   \(e/(a+b)\).

Here \(e\) is the field message whose BCH encoding is all ones. The exact
scan covers all

\[
\binom{16385}{2}=134{,}225{,}920
\]

supports. Both BCH-weight histograms begin at weight 32. Neither histogram
contains weight 22 or weight 128. A direct 12-data-block enumeration matches
the optimized scanner.

The first normalization has 64 fixed `1111` packets and 32 variable
packets. The second has 32 fixed packets and 64 variable packets. Their
respective transfers are

\[
2^{-70.5135}
\quad\text{and}\quad
2^{-55.4175}.
\]

The two normalization families are disjoint in the full scan. Their sum is

\[
\mathbb E Z_D^{p_1,\mathrm{end}}
\le 2^{-55.4174}.
\tag{3}
\]

## Closed subcase and limitation

Equations (2) and (3) cover every outer word of symbol weight three that
contains an all-one BCH block. Together they give

\[
\boxed{
\mathbb E Z_D^{\mathrm{weight\ three,all\ one}}
\le 2^{-44.5519}.
}
\tag{4}
\]

Equation (4) closes the complete all-one endpoint family at distance
\(0.09L\). It does not cover weight-three outer words without an all-one BCH
block. It also does not cover outer occupation four or larger.

## Reproduction

```powershell
python scripts/analyze_riffle_shiftalpha64_endpoint_fixed_packets.py --optimizer-maxiter 100
rustc -O -C target-cpu=native scripts/scan_riffle_shiftalpha64_p1_endpoints.rs -o scripts/scan_riffle_shiftalpha64_p1_endpoints.exe
python scripts/audit_riffle_shiftalpha64_p1_endpoints.py --data-blocks 12
scripts/scan_riffle_shiftalpha64_p1_endpoints.exe 16384
python scripts/analyze_riffle_shiftalpha64_p1_endpoint.py --optimizer-maxiter 120
```

The new artifact hashes are:

- fixed-packet transfer: `46B9F9FF1AE2CBEA1E6BFACAEFBA6CD9C6B1117C8E0F3E68A32743B8AD84A3F9`;
- parity transfer: `65235F7000DA5140113560E3EE164359356802955B4F69B9532C43EDFC9397AD`;
- scanner source: `9FC5BD4F63D38B88255BD6A4B06241A0E0E51FE947C6E7502053949E4FB565B6`;
- scanner executable: `0AE2F2D16CDC15FB9478EC55C7CC3587593748170E2F180C3153E78E09D87AAE`;
- direct audit: `E2730C7B1477AB74371EABE9D18BBE8F07F2380AD89F919E54979EB2949C80A2`.

## Next goal

Remove the all-one endpoint family from the occupation-three sum. Recompute
the remaining dominant weight pattern before designing another compatibility
table. The next calculation should decide whether the remaining
occupation-three shell needs one more exact small-support argument or already
admits a bulk spectral bound.
