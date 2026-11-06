"""Exact rational checks for the sparse Golay-source packet-tail proof."""

from fractions import Fraction


def main() -> None:
    rho = Fraction(4, 25)
    eta = Fraction(1, 5000)
    c0 = 4 * rho / (1 - rho)
    assert c0 == Fraction(16, 21)

    # zeta=sqrt(16/21)<7/8, hence kappa=4*zeta/(1-zeta)<28.
    zeta_upper = Fraction(7, 8)
    assert c0 < zeta_upper**2
    kappa_upper = 4 * zeta_upper / (1 - zeta_upper)
    assert kappa_upper == 28

    # e<68/25. This rational upper bound certifies C_* eta^3<0.59.
    e_upper = Fraction(68, 25)
    ceta3_upper = (Fraction(4095, 24) * e_upper
                    * (4 * kappa_upper / (1 - 8 * eta))**4 * eta**3)
    assert ceta3_upper < Fraction(59, 100)

    # At x=12 eta, y=kappa*x/(1-2x)<e^-2, so ln(y)<-2.
    # The remaining derivative terms satisfy 1+2x/(1-2x)<2.
    x = 12 * eta
    y_upper = kappa_upper * x / (1 - 2*x)
    assert y_upper < 1 / e_upper**2
    assert 2*x/(1-2*x) < 1
    print({"rho": str(rho), "eta": str(eta),
           "ceta3_rational_upper": str(ceta3_upper),
           "sparse_tail": "O(b^-3/2) from paper's accumulator inequality"})


if __name__ == "__main__":
    main()
