from XBasis import *

def singleValuedKZ3(a,b,c,d,mpl=None):
    res = KZCoefDepth3(a,b,c,d) * pi**(a+b+c+d+3)
    res = simplifyWithMaple(mpl,res)
    res = mapToSingleValued(mpl,res)
    return res

def tauDepth3(a,b,c,d):
    result = 0
    for gamma in range(int(m / 3 + 1), m + 1):
        for beta in range(min(m - gamma + 1, gamma)-1,max(0,m-2*gamma)-1,-1):
            alpha = m - beta - gamma
            assert 0 <= beta < gamma and 0 <= alpha <= gamma
            cof = binom(alpha,a) * binom(beta,a+b-alpha) * binom(gamma,d) * Integer(-1)**(gamma-alpha)
            cof -= binom(alpha,d) * binom(beta,a) * binom(gamma,a+b-beta) * Integer(-1)**(alpha-beta)
            cof -= binom(alpha,a) * binom(beta,d) * binom(gamma,a+b-alpha) * Integer(-1)**(beta-alpha)
            cof += binom(alpha,d) * binom(beta,a+b-gamma) * binom(gamma,a) * Integer(-1)**(alpha-gamma)
            result += cabc(alpha,beta,gamma) * Integer(-1)**(2*a+b-d) * cof

    return result

def singleValuedKZ2(a,b,c,mpl=None):
    res = KZCoefDepth2(a,b,c) * pi**(a+b+c+2)
    res = simplifyWithMaple(mpl,res)
    res = mapToSingleValued(mpl,res)
    return res

def tauDepth2(a, b, c):
    result = 0
    n = int((a+b+c+1)/2)
    for beta in range(n, 2*n):
        alpha = 2*n-1-beta
        assert 0 <= alpha < beta
        cof = binom(alpha,a) * binom(beta,c) * Integer(-1)**(a+beta-c)
        cof -= binom(alpha,c) * binom(beta,a) * Integer(-1)**(a+alpha-c)
        result += cof * cab(n,alpha, beta)

    return result
if __name__ == "__main__":
    n = 5
    m = 2*n-1
    mpl = startHyperlogProd()

    for a in range(0, m + 1):
        for b in range(0, m + 1 - a):
            #for d in range(0, m + 1 - a - b):
            c = m - a - b
            assert c >= 0 and a + b + c == m

            svKZ = singleValuedKZ2(a,b,c,mpl)
            tau = tauDepth2(a,b,c)
            tau = simplifyWithMaple(mpl,tau)

            if tau != 0:
                print(svKZ/tau)
            else:
                print("0", svKZ)
