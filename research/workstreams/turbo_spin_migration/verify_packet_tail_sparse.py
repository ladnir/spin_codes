"""Exact rational checks for the sparse Golay-source packet-tail proof."""

from fractions import Fraction


def main() -> None:
    e_upper = Fraction(68, 25)
    profiles = [
        (Fraction(4, 25), Fraction(1, 5000), Fraction(7, 8), 28),
        (Fraction(17, 100), Fraction(1, 10000), Fraction(91, 100), 41),
        (Fraction(9, 50), Fraction(1, 20000), Fraction(15, 16), 60),
        (Fraction(19, 100), Fraction(1, 50000), Fraction(31, 32), 124),
    ]
    for rho, eta, zeta_upper, kappa_upper in profiles:
        c0 = 4 * rho / (1 - rho)
        assert c0 < zeta_upper**2
        assert 4*zeta_upper/(1-zeta_upper) <= kappa_upper
        ceta3_upper = (Fraction(4095, 24) * e_upper
                        * (4 * kappa_upper / (1 - 8 * eta))**4 * eta**3)
        assert ceta3_upper < 1

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
