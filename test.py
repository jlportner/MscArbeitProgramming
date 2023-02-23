import numpy as np
import math
from singleValuedKZ import *

def binom(n,k):
    if k < 0 or n < 0:
        return 0
    else:
        return math.comb(n,k)

"""
#n = 10
mpl = startHyperlogProd()


for n in range(4,5):
    #attention a is not 0
    for a in range(1, 2 * n):
        for c in range(0, 2 * n - a):
            b = 2 * n - 1 - c - a

            p1 = 0
            #for i in range(0,n):
            #    p1 += (-1)**(a-c) * (binom(2*n,i+1) - binom(2*n,i) + (-1)**i) * (-1)**(i+1) * (binom(i,a) * binom(2*n-1-i,c) + binom(i,c) * binom(2*n-1-i,a))

            #for i in range(0,b+1):
                #p1 += (-1)**(a-c) * (binom(2*n,a+i+1) - binom(2*n,a+i) + (-1)**(a+i)) * (-1)**(a+i+1) * (binom(a+i,a) * binom(b+c-i,c))
            #    p1 += (-1)**(a-c) * binom(2*n,a+i+1) * (-1)**(a+i+1) * (binom(a+i,a) * binom(b+c-i,c))

            p1 += (-1)**(a+1) * binom(2*n,a) + (-1)**b * binom(2*n,b) + (-1)**(c+1) * binom(2*n,c)

            p2 = 0
            for i in range(0,a+c+1):
                #p2 += binom(a+i,a) * binom(b+c-i,b) * (-1)**(b-i+1) * (binom(2*n,b+c-i) + binom(2*n,a+i) + (-1)**(b+c-i+1)) * (-1/2)
                #p2 += binom(a+i,a) * binom(b+c-i,b) * (-1)**(b-i+1) * (binom(2*n+1,b+c-i+1) + (-1)**(b+c-i+1)) * (-1/2)
                #p2 += (-1)**c / 2 * binom(a+c-i,a) * binom(b+i,i) * (binom(2*n+1,b+i+1) * (-1)**(b+i) - 1)
                p2 += (-1)**(b+c+i) / 2 * binom(a+c-i,a) * binom(b+i,b) * binom(2*n+1,b+i+1)

            p2 -= (-1)**c /2 * binom(2*n,c)
            if p2 != 0:
                print(p1/p2)
            else:
                print(p1,p2)

"""
def ATDepth2(a,b,c):
    if (a+b+c) % 2 == 1:
        return 0
    n = int((a+b+c+2) / 2)
    result = KZCoefDepth2(a,b,c)
    for l in range(1,n-1):
        m = n-1-l
        interCoef = binom(2*l,a) * binom(2*m,c) * Integer(-1)**(a-c)
        interCoef += binom(2 * l, a) * binom(2 * m, b) * Integer(-1) **(a+b)
        interCoef -= binom(2 * l, c) * binom(2 * m, b) * Integer(-1) **(b+c)
        result += J2lm(l,m,Rational(1/2)) * c2n(l) * c2n(m) * interCoef
    for s in range(1,b+a):
        p = b + a - 2 * s
        interCoef = binom(p+c,p) * binom(2*s,a) * Integer(-1)**(a+c+1)
        interCoef -= binom(p+c,p) * binom(2*s,b) * Integer(-1)**(b+c+1)
        result += zeta(2*s+1) * zeta(p+c+1) / (2*pi*i)**(2*n) * interCoef

    for s in range(1,b+c):
        q = b + c - 2 * s
        interCoef = binom(q+a,q) * binom(2*s,b) * Integer(-1)**(b+q+1)
        result += zeta(2*s+1) * zeta(a+q+1) / (2*pi*i)**(2*n) * interCoef

    return result

mpl = startHyperlogProd()

res = ATDepth2(0,4,2)
print(simplifyWithMaple(mpl,res))