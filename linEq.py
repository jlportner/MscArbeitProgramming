import numpy as np
import math
from mzvCalc import *
from timeit import default_timer as timer

def matLine(n,a,b,c):
    line = []
    for beta in range(n, 2*n):
        alpha = 2*n-1-beta
        assert 0 <= alpha < beta
        cof = binom(alpha,a) * binom(beta,c) * (-1)**(a+beta-c)
        cof -= binom(alpha,c) * binom(beta,a) * (-1)**(a+alpha-c)
        line += [cof]
    return line

def rhsLine(n,a,b,c):
    result = 0

    for s in range(1,int((a+b)/2)+1):
        p = a+b - 2*s
        assert p >= 0
        result += I1(s) * c2n(s) * (binom(2*s,a) * Integer(-1)**a - binom(2*s,b) * Integer(-1)**b) * KZCoef(p,c)

    for s in range(1, int((b+c) / 2) + 1):
        q = b+c - 2 * s
        assert q >= 0
        result += I1(s) * c2n(s) * binom(2*s,b) * Integer(-1)**b * KZCoef(a, q)


    return 1/ I1(n) * (-2 * KZCoef(a,b,c) - result)

def fullMat(n,mpl=None):
    if mpl is None:
        mpl = startHyperlogProd()
    M = []
    rhs = []
    for a in range(0,2*n):
        for c in range(0,2*n-a):
            b = 2*n - 1 - c - a
            mL = matLine(n, a, b, c)
            rhsL = rhsLine(n, a, b, c)
            rhsL = simplifyWithMaple(mpl, rhsL)
            M += [mL]
            rhs += [rhsL]
    M = np.array(M,dtype=int)
    return M,rhs

def fullRankReduction(n,mpl=None):
    if mpl is None:
        mpl = startHyperlogProd()
    M = []
    rhs = []
    varList = []
    b=0
    for c in range(n, 2 * n):
        a = 2 * n - 1 - c
        mL = matLine(n, a, b, c)
        rhsL = rhsLine(n,a,b,c)
        rhsL = simplifyWithMaple(mpl, rhsL)
        M += [mL]
        rhs += [rhsL]
        varList += ["c_{%d,%d}" % (b, c)]

    M = np.array(M,dtype=int)
    return M,rhs,varList

def cDepth2(n,mpl=None):
    M,rhs,varList = fullRankReduction(n,mpl=mpl)
    A = matrix(M)
    b = matrix(rhs).T
    return A.solve_right(b),varList

mpl = startHyperlogProd()
#for n in range(2,10):
#    M,rhs = fullMat(n,mpl)
#    A = matrix(M)
#    b = matrix(rhs).T
#    print("### %d ###" % n)
#    print(A.solve_right(b))

for n in range(2,10):
    print("### %d ###" % n)
    for beta in range(n, 2*n):
        alpha = 2*n-1-beta
        val = cab(n,alpha,beta)
        val = simplifyWithMaple(mpl,val)
        print(val)

#sol,varList = cDepth2(n,mpl=mpl)
#for i,x in enumerate(sol):
#    print(varList[i],x)


