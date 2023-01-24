import warnings
from sage.all import *
#from mzvCalc import *


zeta = function("zeta")
zetaSV = function("zetaSV", latex_name="\\zeta_{\\text{sv}}")
def binom(n,k):
    if k < 0:
        return 0
    else:
        return binomial(n,k)

def evenZeta(n):
    return (-1)**(int(n/2)+1) * (2*pi)**n / (2*factorial(n)) * bernoulli(n)

def parityThm(n,m):
    k = m+n
    if k % 2 != 1:
        #warnings.warn("Parity Theorem only works for odd n+m, returning originial zeta")
        return zeta(n,m)
    else:
        K = int((m+n-1)/2)

    result = 0
    for s in range(0,K):
        result += (binom(k-2*s-1,m-1) + binom(k-2*s-1,n-1) - kronecker_delta(n,2*s) + (-1)**m * kronecker_delta(s,0)) * evenZeta(2*s) * zeta(k-2*s)
        break

    return (-1)**m * result

def zetaDepth1(a,b):
    return (-1)**b * binom(a+b,a) * zeta(a+b+1)

def zetaDepth2(a,b,c):
    result = 0
    delta=0
    if a == 0:
        delta = 1
    for k in range(delta,c+1):
        result += binom(a+k,a) * binom(b+c-k,b) * parityThm(a+k+1,b+c-k+1)

    if a==0:
        accum =0
        for i in range(1,b+c):
            accum += zeta(i+1,b+c-i+1)
        accum += 2*zeta(b+c+1,1)
        result -= binom(b+c,b) * accum
    return (-1)**c * result

def zetaDepth3(a,b,c,d):
    result = 0
    delta=0
    if a == 0:
        delta = 1
    for k in range(delta,d + 1):
        for j in range(d-k+1):
            result += binom(a+k,a) * binom(b + j, b) * binom(c+d - k-j, c) * zeta(a+k+1,b+j+1,c+d-j-k+1)

    if a==0:
        delta = 0
        if b==0:
            delta = 1
        for j in range(delta,d+1):
            accum= 0
            for i in range(1,b+j+1):
                accum -= zeta(i+1,b+j-i+1,c+d-j+1)
            for i in range(0,c+d-j+1):
                accum -= zeta(b+j+1,i+1,c+d-j-i+1)

            accum -= zeta(b+j+1,c+d-j+1,1)
            result += binom(b + j, b) * binom(c+d-j, c) * accum

        if b==0:
            accum = 0
            for i in range(1,c+d):
                for j in range(1,i+1):
                    accum += zeta(j+1,i-j+1,c+d-i+1)
                for j in range(0,c+d-i):
                    accum += zeta(i+1,j+1,c+d-i-j+1)
                accum += 4* zeta(i+1,c+d-i+1,1)
            accum += 6* zeta(c+d+1,1,1)

            result += Rational(1/2) * binom(c+d,c) * accum
    return (-1)**d * result

def convertToBrownBasis(expr):
    z533to353 = Rational(1/378)*pi**6*zeta(5) - Rational(1/30)*pi**4*zeta(7) - Rational(15/2)*pi**2*zeta(9) + Rational(1/2)*zeta(5, 3)*zeta(3) + Rational(299/4)*zeta(11) - Rational(1/2)*zeta(3,5,3)
    z553to535 = Rational(1/3780)*pi**8*zeta(5) - Rational(5/18)*pi**4*zeta(9) - Rational(275/12)*pi**2*zeta(11) + Rational(1/2)*zeta(5, 3)*zeta(5) + Rational(1003/4)*zeta(13) - Rational(1/2) * zeta(5,3,5)
    zeta733to373 = Rational(1/3150)*pi**8*zeta(5) + Rational(4/945)*pi**6*zeta(7) - Rational(14/45)*pi**4*zeta(9) - Rational(407/12)*pi**2*zeta(11) + Rational(1/2)*zeta(7, 3)*zeta(3) + Rational(358)*zeta(13) - Rational(1/2) * zeta(3,7,3)

    expr = expr.subs(zeta(5,3,3)==z533to353)
    expr = expr.subs(zeta(5,5,3)==z553to535)
    expr = expr.subs(zeta(7,3,3)==zeta733to373)
    return expand(expr)

def convertToSingleValued(expr):
    z353tozSV353 = Rational(1/2) * zetaSV(3,5,3) + zeta(3) * zeta(5,3) + 5 * zeta(3)**2 * zeta(5)
    z535tozSV535 =  Rational(1/2) * zetaSV(5,3,5) + 11*zeta(5)*zeta(5,3) + 60*zeta(5)**2 * zeta(3) + 5*zeta(5) * evenZeta(8)
    z373tozSV373 = Rational(1/2) * zetaSV(3,7,3) + zeta(3) * zeta(7,3) + 14 * zeta(3)**2 * zeta(7) + 12*zeta(5) * zeta(5,3) + 72 * zeta(5)**2 * zeta(3) + 6*zeta(5)*evenZeta(8)
    expr = expr.subs(zeta(3,5,3)==z353tozSV353)
    expr = expr.subs(zeta(5,3,5)==z535tozSV535)
    expr = expr.subs(zeta(3,7,3)==z373tozSV373)
    return expand(expr)

def convertToSingleValuedSchnetz(expr):
    z533tozSV533 = Rational(1/2) * (Rational(1/189)*pi**6*zeta(5) - Rational(1/15)*pi**4*zeta(7) - 15*pi**2*zeta(9) + 5*zeta(5)*zeta(3)**2 + zetaSV(5, 3, 3))
    z553tozSV553 = Rational(1/2) * (Rational(-5/9)*pi**4*zeta(9) - Rational(275/6)*pi**2*zeta(11) - 50*zeta(5)**2*zeta(3) - 10*zeta(5, 3)*zeta(5) + zetaSV(5, 5, 3))
    z733tozSV733 = Rational(1/2) * (Rational(8/945)*pi**6*zeta(7) - Rational(28/45)*pi**4*zeta(9) - Rational(407/6)*pi**2*zeta(11) - 60*zeta(5)**2*zeta(3) + 14*zeta(7)*zeta(3)**2 - 12*zeta(5,3)*zeta(5) + zetaSV(7,3,3))
    z933tozSV933 = Rational(1/2) * (Rational(-1/1575)*pi**8*zeta(7) - Rational(29/945)*pi**6*zeta(9) - Rational(14/5)*pi**4*zeta(11) - Rational(403/2)*pi**2*zeta(13) - 72*zeta(5)**3 - 318*zeta(7)*zeta(5)*zeta(3) + 27*zeta(9)*zeta(3)**2 - 30*zeta(7)*zeta(5, 3) - 12*zeta(7, 3)*zeta(5) + zetaSV(9,3,3))
    z735tozSV735 = Rational(1/2) * (Rational(2/31185)*pi**10*zeta(5) + Rational(2/675)*pi**8*zeta(7) + Rational(34/945)*pi**6*zeta(9) - Rational(11/9)*pi**4*zeta(11) - Rational(1001/6)*pi**2*zeta(13) + 78*zeta(5)**3 + 336*zeta(7)*zeta(5)*zeta(3) + 28*zeta(7)*zeta(5, 3) + 12*zeta(7, 3)*zeta(5) + zetaSV(7,3,5))
    expr = expr.subs(zeta(5,3,3)==z533tozSV533)
    expr = expr.subs(zeta(5,5,3)==z553tozSV553)
    expr = expr.subs(zeta(7,3,3)==z733tozSV733)
    expr = expr.subs(zeta(9,3,3)==z933tozSV933)
    expr = expr.subs(zeta(7,3,5)==z735tozSV735)
    return expand(expr)

if __name__ == "__main__":
    mpl = startHyperlogProd()

    n = 4
    m = 2*n-2
    #print(zetaDepth2(2,2,3))
    #exit()
    print("Start")
    for a in range(0, m + 1):
        for b in range(0, m + 1 - a):
            for d in range(0, m + 1 - a - b):
                c = m - a - b - d
                #a,b,c,d = 0,0,3,1
                #print(a,b,c,d)
                D1 =(simplifyWithMaple(mpl,regZetaWithHyperlog(a,b,c,d,mpl=mpl)))
                D2 = (simplifyWithMaple(mpl,zetaDepth3(a,b,c,d)))
                D3 = (simplifyWithMaple(mpl, zetaDepth3Correct(a,b,c,d)))
                if D1 != D2 or D1 != D3:
                    print(a,b,c,d)
                #exit()
