import SpinCodes.Structured.MapSpectrumDefs

namespace Spin.Structured.MapSpectrum

def sumHistograms (hs : List (List Nat)) : List Nat :=
  hs.foldr (List.zipWith (· + ·)) (List.replicate 129 0)

end Spin.Structured.MapSpectrum
