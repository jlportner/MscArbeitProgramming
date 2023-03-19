from depth3Formula import *

def adjustToBasis(alpha, beta, gamma, n):
    res = cabc(alpha, beta, gamma)
    res = convertToSingleValuedSchnetz(res)
    if n == 4:
        e1 = var("e1")
        res = res.subs(zeta(9) == e1 * (2 * pi * i) ** 9)
        res = res.expand()
    elif n == 5:
        e1, e2 = var("e1,e2")
        res = res * (2 * pi * I) ** 11
        res = res.subs(
            zetaSV(5, 3, 3) == (Rational(5 / 7759752) * e2 - Rational(22020 / 3553) * zeta(3) ** 2 * zeta(5)),
            zeta(11) == Rational(1 / 116396280) * e1)
        res = res.expand()
        # res = res.subs(e2 == 0)
    elif n == 6:
        e1, e2, e3 = var("f1,f2,f3")
        res = res * (2 * pi * I) ** 13
        res = res.subs(zeta(13) == Rational(1 / 2974571600) * e1,
                       zetaSV(5, 5, 3) == Rational(1 / 13520780) * e2 + Rational(203950 / 5681) * zeta(
                           5) ** 2 * zeta(3))
        res = res.subs(
            zetaSV(7, 3, 3) == Rational(1 / 19315400) * e3 + Rational(244740 / 5681) * zeta(5) ** 2 * zeta(
                3) - Rational(
                123508 / 7429) * zeta(7) * zeta(3) ** 2)
        res = res.expand()
        # res = res.subs(e2 == 0, e3 == 0)
    elif n == 7:
        e1, e2, e3 = var("g1,g2,g3")
        res = res * (2 * pi * I) ** 15
        res = res.subs(zeta(15) == Rational(1 / 70578471600) * e1,
                       zetaSV(9, 3, 3) == Rational(1 / 258529200) * e2 + 48 * zeta(5) ** 3 + Rational(
                           8694314 / 37145) * zeta(7) * zeta(5) * zeta(3) - Rational(68094 / 2185) * zeta(9) * zeta(
                           3) ** 2)
        res = res.subs(
            zetaSV(7, 3, 5) == Rational(7 / 1706292720) * e3 - 56 * zeta(5) ** 3 - Rational(5826772 / 22287) * zeta(
                7) * zeta(5) * zeta(3))
        res = res.expand()
    return res

mpl = startHyperlogProc()
for n in range(3,8):
    print("\\subsection{n=%d}" % n, )
    print("\\begin{align*}")
    m = 2*n-2
    for gamma in range(int(m / 3 + 1), m + 1):
        for beta in range(min(m - gamma + 1, gamma)-1,max(0,m-2*gamma)-1,-1):
            alpha = m - beta - gamma
            assert 0 <= beta < gamma and 0 <= alpha <= gamma
            x = adjustToBasis(alpha,beta,gamma,n)
            latexStr = latex(x)
            print("c_{%d,%d,%d}" % (alpha,beta,gamma), "&=", latexStr, "\\\\")
    print("\\end{align*}")