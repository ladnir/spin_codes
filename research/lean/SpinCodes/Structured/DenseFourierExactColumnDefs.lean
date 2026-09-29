import SpinCodes.Structured.DenseFourierExactDefs

namespace Spin.Structured.DenseFourierExact

structure ICoords where
  Z : Int
  D : Int
  S : Fin 5 → Int

def weightedTotal (v : ICoords) : Int :=
  v.Z+((List.finRange 5).map fun (i : Fin 5) => (SparsePolynomial.counts.getD i 1:Int)*v.S i).sum

def targetD (qn qd zn zd : Int) : Int :=
  max (target qn qd zn zd 48) (max (target qn qd zn zd 56)
    (max (target qn qd zn zd 64) (max (target qn qd zn zd 72) (target qn qd zn zd 80))))

def zeroTerm (qn qd zn zd : Int) (v : ICoords) (j : Nat) : Int :=
  let c : Int := polyChoose 128 j
  let k : Int := (SparsePolynomial.Data.weights.getD j SparsePolynomial.Data.weight0).kernel
  (k*v.Z+(c-k)*v.D)*(qn*zn)^j*((qd-qn)*zd)^(128-j)

def zeroColumn (qn qd zn zd : Int) (v : ICoords) : Int :=
  ((List.range 129).map (zeroTerm qn qd zn zd v)).sum

def checkZ (qn qd zn zd rn rd : Int) (v : ICoords) : Bool :=
  decide (zeroColumn qn qd zn zd v*rd ≤ rn*(qd*zd)^128*v.Z)

def checkD (qn qd zn zd rn rd : Int) (v : ICoords) : Bool :=
  decide (targetD qn qd zn zd*weightedTotal v*rd ≤ rn*denominator qd zd*v.D)

def checkS (qn qd zn zd rn rd : Int) (v : ICoords) (i : Fin 5) : Bool :=
  decide (target qn qd zn zd (SparsePolynomial.levels.getD i 0)*weightedTotal v*rd ≤
    rn*denominator qd zd*v.S i)

end Spin.Structured.DenseFourierExact

