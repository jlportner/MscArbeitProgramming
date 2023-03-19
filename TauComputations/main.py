from ATDepth3 import *

print("Enter coefficient weight: ")
n = int(input())

mpl = startHyperlogProc()
print("c_{%d} = " %(2*n), latex(c2n(n)))

m = 2 * n - 1
for beta in range(m, n - 1, -1):
    alpha = m - beta
    assert 0 <= alpha < beta
    expr = "c_{%d,%d} = " % (alpha, beta)
    print(expr,latex(cab(n,alpha, beta)))

m = 2 * n - 2
for gamma in reversed(range(int(m / 3 + 1), m + 1)):
    for beta in range(min(m - gamma + 1, gamma) - 1, max(0, m - 2 * gamma) - 1, -1):
        alpha = m - beta - gamma
        assert 0 <= beta < gamma and 0 <= alpha <= gamma
        expr = "c_{%d,%d,%d} = " % (alpha, beta, gamma)
        print(expr, latex(cabc(alpha,beta,gamma)))


print("Enter word length for AT-Associator: ")
m = int(input()) - 3
for a in range(0, m + 1):
    for b in range(0, m + 1 - a):
        for d in range(0, m + 1 - a - b):
            c = m - a - b - d
            x = interpolationAssociatorDepth3(a, b, c, d, t=Rational(1 / 2), mpl=mpl)
            latexStr = latex(x)
            word = ""
            for i in (a, b, c, d):
                if i == 0:
                    word += "y "
                else:
                    word += "x^%d y " % i
            word = word[:-2].replace("y y y", "y^3").replace("y y", "y^2")
            print("\Phi_{AT}(%s) = " % word, latexStr)