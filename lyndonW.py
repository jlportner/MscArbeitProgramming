from sage.all import *
import sage.combinat.words.lyndon_word as ly
import numpy as np

L= LieAlgebra(QQ, "x,y", representation="polynomial")
#R=PolynomialRing(QQ,"w,z",commutative=False)
x,y = L.gens()

expansions = []
for w in ly.LyndonWords([3,6]):
    bra = ly.standard_bracketing(w)

    braStr = str(bra).replace("2","x").replace("1","y")
    brak = eval(braStr)
    #bra = SR(braStr)
    print(str(w).replace("2","x").replace("1","y"), L(brak))
    #a = R(str(L(brak)).replace("x","w").replace("y","z"))
    #expansions += [a]

M = []
n=4
m = 2*n-2
for a in range(0,m+1):
    for b in range(0,m+1-a):
        for d in range(0,m+1-a-b):
            c = m - a - b - d
            line = []
            for terms in expansions:
                cof  = terms.coefficients(x**a * y * x**b * y * x**c * y * x**d)
                line += [cof]
            M += [line]

M = np.array(M, dtype=int)
print(M)

def ad(s):
    if s == 0:
        return y
    else:
        expr = [x,y]
        for i in range(s-1):
            expr = [x,expr]

    return expr

#print(ad(5))