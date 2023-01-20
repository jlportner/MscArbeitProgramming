import pickle

import numpy as np

from mzvCalc import *
import sage.misc.persist

def matLine(m,a,b,c,d):
    line = []
    for gamma in range(int(m / 3 + 1), m + 1):
        for beta in range(min(m - gamma + 1, gamma)-1,max(0,m-2*gamma)-1,-1):
            alpha = m - beta - gamma
            assert 0 <= beta < gamma and 0 <= alpha <= gamma
            cof = binom(alpha,a) * binom(beta,a+b-alpha) * binom(gamma,d) * Integer(-1)**(gamma-alpha)
            cof -= binom(alpha,d) * binom(beta,a) * binom(gamma,a+b-beta) * Integer(-1)**(alpha-beta)
            cof -= binom(alpha,a) * binom(beta,d) * binom(gamma,a+b-alpha) * Integer(-1)**(beta-alpha)
            cof += binom(alpha,d) * binom(beta,a+b-gamma) * binom(gamma,a) * Integer(-1)**(alpha-gamma)
            line += [Integer(-1)**(2*a+b-d) * cof]
    return line

def fullMatAndRhs(n,mpl=None):
    if mpl is None:
        mpl = startHyperlogProd()
    m = 2 * n - 2
    M = []
    rhs = []
    for a in range(0,m+1):
        for b in range(0,m+1-a):
            for d in range(0,m+1-a-b):
                c = m - a - b - d
                assert c >= 0 and a+b+c+d == m
                mL =  matLine(m, a, b, c, d)
                #print(a,b,c,d,mL)
                line = rhsLine(n, a, b, c, d)
                line = simplifyWithMaple(mpl, line)
                M += [mL]
                rhs += [line]
                #print(mL, line)

    M = np.array(M, dtype=int)
    return M, rhs

def fullRankMatAndRhs(n,mpl=None):
    if mpl is None:
        mpl = startHyperlogProd()
    m = 2*n-2
    M = []
    varList = []
    rhs = []
    for d in range(int(m / 3 + 1), m + 1):
        for b in range(min(m - d + 1, d) -1 , max(0,m-2*d)-1,-1):
            a = m - b - d
            c = 0
            varList += [(a,b,d)]
            M += [matLine(m, a, b, c, d)]
            line = simplifyWithMaple(mpl,rhsLine(n,a,b,c,d))
            rhs += [line]


    M = np.array(M, dtype=int)
    return M, rhs, varList

def rhsDoubleIntegralPart(n,a,b,c,d, t=1):
    result = 0
    for l in range(1,int((a+b+c)/2.0)):
        for m in range(1,int((a+b+c)/2.0) - l + 1):
            p = a+b+c-2*l-2*m
            assert p >= 0
            interCoef = binom(2*l,a) * binom(2*m,a+b-2*l) * Integer(-1)**b
            interCoef -= binom(2*l,a) * binom(2*m,c) * Integer(-1)**(a-c)
            interCoef += binom(2*l,b) * binom(2*m,a) * Integer(-1)**(a+b)
            interCoef -= binom(2*l,b) * binom(2*m,a+b-2*l) * Integer(-1)**a
            interCoef -= binom(2*l,c) * binom(2*m,a) * Integer(-1)**(a-c)
            interCoef -= binom(2*l,b) * binom(2*m,b+c-2*l) * Integer(-1)**c
            interCoef += binom(2*l,b) * binom(2*m,c) * Integer(-1)**(b+c)
            interCoef += binom(2*l,c) * binom(2*m,b+c-2*l) * Integer(-1)**b
            result += J2lm(l,m,t) * c2n(l) * c2n(m) * KZCoefDepth1(p,d) * interCoef

    for l in range(1,int((b+c+d)/2.0)):
        for m in range(1,int((b+c+d)/2.0) - l + 1):
            q = b+c+d-2*l-2*m
            assert q >= 0
            interCoef = binom(2*l,b) * binom(2*m,b+c-2*l) * Integer(-1)**c
            interCoef += binom(2*l,c) * binom(2*m,b) * Integer(-1)**(b+c)
            interCoef -= binom(2*l,c) * binom(2*m,b+c-2*l) * Integer(-1)**b
            result += J2lm(l,m,t) * c2n(l) * c2n(m) * KZCoefDepth1(a,q) * interCoef

    for l in range(1,int((c+d)/2.0) + 1):
        for m in range(1,int((a+b)/2.0) + 1):
            p = a+b-2*m
            q = c+d-2*l
            assert p >= 0 and q >= 0
            interCoef = binom(2*l,c) * binom(2*m,a) * Integer(-1)**(a+c)
            interCoef -= binom(2*l,c) * binom(2*m,b) * Integer(-1)**(c-b)
            result += J2lm(l,m,t) * c2n(l) * c2n(m) * KZCoefDepth1(p,q) * interCoef

    for l in range(1,int((a+b)/2.0) + 1):
        for m in range(1,int((c+d)/2.0) + 1):
            p = a+b-2*l
            q = c+d-2*m
            assert p >= 0 and q >= 0
            interCoef = binom(2*l,a) * binom(2*m,c) * Integer(-1)**(a+c)
            interCoef -= binom(2*l,b) * binom(2*m,c) * Integer(-1)**(c-b)
            result += J2lm(l,m,t) * c2n(l) * c2n(m) * KZCoefDepth1(p,q) * interCoef

    return result

def rhsTripleIntegralPart(n,a,b,c,d, t=1):
    n = int((a+b+c+d +2)/2)
    result = 0
    for l in range(1,(n-3)+1):
        for m in range(1,(n-2-l)+1):
            h = (n-1)-m-l
            assert h >= 1
            cof = binom(2*l,a) * binom(2*m,a+b-2*l) * binom(2*h,d) * Integer(-1)**(b-d)
            cof +=binom(2*l,a) * binom(2*m,c) * binom(2*h,a+b-2*l) * Integer(-1)**(b+c)
            cof -=binom(2*l,a) * binom(2*m,c) * binom(2*h,d) * Integer(-1)**(a-c-d)

            cof +=binom(2*l,b) * binom(2*m,a) * binom(2*h,d) * Integer(-1)**(a+b-d)
            cof -=binom(2*l,b) * binom(2*m,a+b-2*l) * binom(2*h,d) * Integer(-1)**(a-d)
            cof +=binom(2*l,c) * binom(2*m,a) * binom(2*h,a+b-2*m) * Integer(-1)**(b+c)
            cof -=binom(2*l,c) * binom(2*m,a) * binom(2*h,d) * Integer(-1)**(a-c-d)
            cof +=binom(2*l,b) * binom(2*m,b+c-2*l) * binom(2*h,a) * Integer(-1)**(a+c)
            cof -=binom(2*l,b) * binom(2*m,c) * binom(2*h,a+b-2*l) * Integer(-1)**(a+c)
            cof +=binom(2*l,c) * binom(2*m,b) * binom(2*h,a) * Integer(-1)**(a+b+c)
            cof -=binom(2*l,c) * binom(2*m,b+c-2*l) * binom(2*h,a) * Integer(-1)**(a+b)
            cof -=binom(2*l,b) * binom(2*m,a+d-2*h) * binom(2*h,d) * Integer(-1)**(a+b)
            cof +=binom(2*l,b) * binom(2*m,c) * binom(2*h,d) * Integer(-1)**(b+c+d)
            cof -=binom(2*l,c) * binom(2*m,b) * binom(2*h,a+b-2*m) * Integer(-1)**(a+c)
            cof +=binom(2*l,c) * binom(2*m,a+d-2*h) * binom(2*h,d) * Integer(-1)**(a-c)

            result += cof * K3hlm(l,m,h,t) * c2n(h) * c2n(l) * c2n(m)

    return result

def rhsCabPart(n,a,b,c,d, t=1):
    result = 0
    for s in range(1, int((a + b + c + 1) / 2) + 1):
        p = a + b + c - (2*s-1)
        assert p >= 0
        for beta in range(s, 2 * s):
            alpha = 2 * s - 1 - beta
            assert 0 <= alpha < beta
            interCoef = binom(alpha,a) * binom(beta,a+b-alpha) * Integer(-1)**(b-alpha)
            interCoef -= binom(alpha,a+b-beta) * binom(beta,a) * Integer(-1)**(b-beta)

            interCoef -= binom(alpha,b+c-beta) * binom(beta,c) * Integer(-1)**(b-alpha)
            interCoef += binom(alpha,c) * binom(beta,b+c-alpha) * Integer(-1)**(b-beta)
            result += I1(s,t) * cab(s,alpha,beta) * interCoef * KZCoefDepth1(p,d) * (-1)

    for s in range(1, int((b + c + d + 1) / 2) + 1):
        q = b + c + d - (2*s-1)
        assert q >= 0
        for beta in range(s, 2 * s):
            alpha = 2 * s - 1 - beta
            assert 0 <= alpha < beta
            interCoef = binom(alpha,b) * binom(beta,b+c-alpha) * Integer(-1)**(c-alpha)
            interCoef -= binom(alpha,b+c-beta) * binom(beta,b) * Integer(-1)**(c-beta)
            result += I1(s,t) * cab(s,alpha,beta) * interCoef * KZCoefDepth1(a,q) * (-1)

    return result

def rhsCDoubleKZpart(n,a,b,c,d,t=1):
    result = 0
    for s in range(1, int((a+b) / 2) + 1):
        p = a+b - 2 * s
        assert p >= 0
        interCoef = binom(2*s,a) * Integer(-1)**a
        interCoef -= binom(2*s,b) * Integer(-1)**b
        result += I1(s,t) * c2n(s) * interCoef * KZCoefDepth2(p,c,d)

    for s in range(1, int((b+c) / 2) + 1):
        q = b+c - 2 * s
        assert q >= 0
        interCoef = binom(2*s,b) * Integer(-1)**b
        interCoef -= binom(2*s,c) * Integer(-1)**c
        result += I1(s,t) * c2n(s) * interCoef * KZCoefDepth2(a,q,d)

    for s in range(1, int((c + d) / 2) + 1):
        r = c + d - 2 * s
        assert r >= 0
        interCoef = binom(2*s,c) * (-1)**c
        result += I1(s,t) * c2n(s) * interCoef * KZCoefDepth2(a,b,r)

    return result

def rhsLine(n,a,b,c,d):
    cabContribution = rhsCabPart(n,a,b,c,d)
    cDoubleKZContribution = rhsCDoubleKZpart(n,a,b,c,d)
    doubleContribution = rhsDoubleIntegralPart(n,a,b,c,d)
    tripleContribution = rhsTripleIntegralPart(n,a,b,c,d)
    return (- 2*KZCoefDepth3(a,b,c,d) - cDoubleKZContribution - cabContribution -doubleContribution - tripleContribution) / I1(n)

def cabcPart(n,a,b,c,d,t=1):
    result = 0
    m = 2*n-2
    for gamma in range(int(m / 3 + 1), m + 1):
        for beta in range(min(m - gamma + 1, gamma)-1,max(0,m-2*gamma)-1,-1):
            alpha = m - beta - gamma
            assert 0 <= beta < gamma and 0 <= alpha <= gamma
            cof = binom(alpha,a) * binom(beta,a+b-alpha) * binom(gamma,d) * Integer(-1)**(gamma-alpha)
            cof -= binom(alpha,d) * binom(beta,a) * binom(gamma,a+b-beta) * Integer(-1)**(alpha-beta)
            cof -= binom(alpha,a) * binom(beta,d) * binom(gamma,a+b-alpha) * Integer(-1)**(beta-alpha)
            cof += binom(alpha,d) * binom(beta,a+b-gamma) * binom(gamma,a) * Integer(-1)**(alpha-gamma)
            result += Integer(-1)**(2*a+b-d) * cof * cabc(alpha,beta,gamma) * I1(n,t)

    return result

def createCache():
    mpl = startHyperlogProd()

    cache = {}
    for n in range(2,7):
        c, varlist = cDepth3(n,mpl)
        for val,key in zip(c,varlist):
            cache[key] = val

    data = str(cache)
    data = data.replace("^","**")
    with open('cabcCache.dat', 'w') as f:
        f.write(data)
    return




def cabc(alpha,beta,gamma):
    if not hasattr(cabc, "cache"):
        with open('cabcCache.dat', 'r') as f:
            data = f.read()
            cabc.cache = sage_eval(data,locals={'zeta': zeta})


    return cabc.cache[(alpha,beta,gamma)]


def cDepth3(n,mpl=None):
    A,b,varList = fullRankMatAndRhs(n,mpl)
    A = matrix(A)
    #print(A)
    b = matrix(b).T
    sol = A.solve_right(b)
    sol = [x[0] for x in sol]
    return sol,varList


def checkDimension(exprList):
    base = [zetaSV(5,5,3), zetaSV(7,3,3), zeta(5)**2 * zeta(3), zeta(7) * zeta(3)**2, zeta(13)]
    #base = [zetaSV(5,3,3),zeta(5) * zeta(3)**2, zeta(11),]
    M = []
    for expr in exprList:
        line = []
        expr = convertToSingleValuedSchnetz(expr)
        expr = SR(str(expr))
        for b in base:
            coef = expr.coefficient(b)
            line += [coef]
            expr -= coef * b
        if expr != 0:
            raise ValueError("Part left over:" + expr)
        M += [line]

    M = matrix(SR,M)
    return M

if __name__ == "__main__":
    #createCache()




    #x^3 y y y x^5
    n = 6
    m = 2 * n - 2
    mpl = startHyperlogProd()
    sol, _ = cDepth3(n,mpl)
    M = checkDimension(sol)
    print(M)
    print("Rank:", M.rank())
    print(M.echelon_form())
    exit()
    expr = simplifyWithMaple(mpl,rhsLine(5,2,0,0,6))
    expr2 = (convertToBrownBasis(expr))
    expr2 = simplifyWithMaple(mpl,expr2)
    print(expr-expr2)
