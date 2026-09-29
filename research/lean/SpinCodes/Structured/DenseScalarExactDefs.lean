namespace Spin.Structured.DenseScalarExact

/-- Exact rational bases, retaining tiny scalar contractions without fixed absolute rounding. -/
def g0 (qn qd zn zd : Int) : Int := (qd-qn)*zd+qn*zn
def g1 (qn qd zn zd : Int) : Int := qn*zd+(qd-qn)*zn

def check (qn qd zn zd rn rd : Int) (d : Nat) : Bool :=
  decide (g0 qn qd zn zd^(128-d)*g1 qn qd zn zd^d*rd ≤ rn*(qd*zd)^128)

end Spin.Structured.DenseScalarExact
