def matLineExtra(m, a, b, c, d):
    line = []
    for gamma in range(int(m / 3 + 1), int(m/2) + 1):
            alpha = gamma
            beta = m - 2*gamma
            cof = binom(alpha,a) * binom(beta,(b+a)-alpha) * binom(gamma,d) * (-1)**(2*a+b-d + gamma - alpha)
            cof2 = binom(alpha,d) * binom(beta,a) * binom(gamma,(b+a)-beta) * (-1)**(2*a+b-d + alpha - beta)
            cof2 -= binom(alpha,a) * binom(beta,d) * binom(gamma,(b+a)-alpha) * (-1)**(2*a+b-d + alpha - beta)
            cof2 += binom(alpha,d) * binom(beta,(b+a)-gamma) * binom(gamma,a) * (-1)**(2*a+b-d + alpha - gamma)
            line += [cof + cof2]
    return line

    """
    for d in range(int(m / 3 + 1), int(m / 2) + 1):
        a,b,c = d,m - 2*d,0
        print(a, b, c, d)
        M += [matLine(m, a, b, c, d) + matLineExtra(m, a, b, c, d)]
    """

def rhsCabPart2(n,a,b,c,d):
    result = 0
    for s in range(1, int((a + b + c + 1) / 2) + 1):
        p = a + b + c - (2*s-1)
        assert p >= 0
        for beta in range(s, 2 * s):
            alpha = 2 * s - 1 - beta
            interCoef = binom(alpha, a) * binom(beta, b + a - alpha) * (Integer(-1)) ** (2 * a + b - alpha)
            interCoef -= binom(alpha, b + a - beta) * binom(beta, a) * (Integer(-1)) ** (2 * a + b - beta)

            interCoef += binom(alpha, b + c - beta) * binom(beta, c) * (Integer(-1)) ** (2 * c + b - beta)
            interCoef -= binom(alpha, c) * binom(beta, b + c - alpha) * (Integer(-1)) ** (2 * c + b - alpha)
            result += cab(s, alpha, beta) * I1(s) * interCoef * KZCoef(p, d) * (-1)**(alpha+beta) #-1 from alpha+beta being odd

    for s in range(1, int((d + b + c + 1) / 2) + 1):
        q = d + b + c - (2*s-1)
        assert q >= 0
        for beta in range(s, 2 * s):
            alpha = 2 * s - 1 - beta
            interCoef = binom(alpha, b) * binom(beta, b + c - alpha) * (Integer(-1)) ** (2 * b + c - alpha)
            interCoef -= binom(alpha, b + c - beta) * binom(beta, b) * (Integer(-1)) ** (2 * b + c - beta)
            result += cab(s, alpha, beta) * I1(s) * interCoef * KZCoef(a, q) * (-1)**(alpha+beta) #-1 from alpha+beta being odd

    return result

def matLine2(m,a,b,c,d):
    line = []
    for gamma in range(int(m / 3 + 1), m + 1):
        for beta in range(min(m - gamma + 1, gamma)-1,max(0,m-2*gamma)-1,-1):
            alpha = m - beta - gamma
            assert 0 <= beta < gamma and 0 <= alpha <= gamma
            cof = binom(alpha,a) * binom(beta,(b+a)-alpha) * binom(gamma,d) * (-1)**(2*a+b-d + gamma - alpha)
            cof -= binom(alpha,d) * binom(beta,a) * binom(gamma,(b+a)-beta) * (-1)**(2*a+b-d + alpha - beta)
            cof -= binom(alpha,a) * binom(beta,d) * binom(gamma,(b+a)-alpha) * (-1)**(2*a+b-d + beta - alpha)
            cof += binom(alpha,d) * binom(beta,(b+a)-gamma) * binom(gamma,a) * (-1)**(2*a+b-d + alpha - gamma)
            line += [cof * (-1)**(alpha+beta+gamma)]
    return line

    result = 0
    for s in range(1, int((b+c) / 2) + 1):
        q = b+c - 2 * s
        assert q >= 0
        result += I1(s) * c2n(s) * (binom(2*s, b) * (Integer(-1))**b - binom(2*s, c) * (Integer(-1))**c) * KZCoef(a, q, d)

    for s in range(1, int((a + b) / 2) + 1):
        p = a + b - 2 * s
        assert p >= 0
        result += I1(s) * c2n(s) * (binom(2*s,a) * (Integer(-1))**a - binom(2*s, b) * (Integer(-1))**b) * KZCoef(p, c, d)

    for s in range(1, int((c + d) / 2) + 1):
        r = c + d - 2 * s
        assert r >= 0
        result += I1(s) * c2n(s) * (binom(2*s, c) * (Integer(-1))**c) * KZCoef(a, b, r)
def rhsDoubleIntegralPart2(n,a,b,c,d):
    result = 0
    for l in range(1,int((a+b+c)/2) + 1):
        for m in range(1, int((a + b + c) / 2 - l) + 1):
            p = a + b + c - 2 * l - 2 * m
            assert p >= 0
            cof = binom(2*l,a) * binom(2*m,a+b-2*l) * Integer(-1)**b
            cof -= binom(2*l,a+b-2*m) * binom(2*m,a) * Integer(-1)**b
            cof += binom(2*l,b) * binom(2*m,a) * Integer(-1)**(a+b)
            cof -= binom(2*l,b) * binom(2*m,a+b-2*l) * Integer(-1)**a
            cof -= binom(2*l,a) * binom(2*m,b) * Integer(-1)**(a+b)
            cof += binom(2*l,a+b-2*m) * binom(2*m,b) * Integer(-1)**a

            cof -= binom(2*l,b+c-2*m) * binom(2*m,c) * Integer(-1)**(b)
            cof += binom(2*l,c) * binom(2*m,b+c-2*l) * Integer(-1)**(b)
            cof -= binom(2*l,b) * binom(2*m,b+c-2*l) * Integer(-1)**(c)
            cof += binom(2*l,b) * binom(2*m,c) * Integer(-1)**(b+c)
            cof += binom(2*l,b+c-2*m) * binom(2*m,b) * Integer(-1)**(c)
            cof -= binom(2*l,c) * binom(2*m,b) * Integer(-1)**(b+c)
            result += cof * J2lm(l,m) * c2n(l) * c2n(m) * KZCoef(p,d)

    for l in range(1, int((b + c + d) / 2) + 1):
        for m in range(1, int((b + c + d) / 2 - l) + 1):
            q = b + c + d - 2 * l - 2 * m
            cof = binom(2*l,b) * binom(2*m,b+c-2*l) * Integer(-1)**c
            cof -= binom(2*l,b+c-2*m) * binom(2*m,b) * Integer(-1)**c
            cof += binom(2*l,c) * binom(2*m,b) * Integer(-1)**(b+c)
            cof -= binom(2*l,c) * binom(2*m,b+c-2*l) * Integer(-1)**b
            cof -= binom(2*l,b) * binom(2*m,c) * Integer(-1)**(b+c)
            cof += binom(2*l,b+c-2*m) * binom(2*m,c) * Integer(-1)**b
            result += cof * J2lm(l, m) * c2n(l) * c2n(m) * KZCoef(a, q)

    return result


def cab(n, alpha, beta=-1):
    assert n > 0
    if beta == -1:
        beta = 2 * n - 1 - alpha
    elif beta + alpha != 2 * n - 1:
        raise ValueError("alpha + beta = 2n-1 needs to be satisfied!")
    if alpha >= beta:
        raise ValueError("alpha < beta needs to be satisfied!")
    a = alpha
    c = beta
    accum = 0

    for s in range(1, int(a / 2) + 1):
        p = a - 2 * s
        assert p >= 0
        accum += I1(s) * c2n(s) * (binom(2 * s, a) * (-1) ** a - 1) * KZCoef(p, c)

    for s in range(1, int(c / 2) + 1):
        q = c - 2 * s
        assert q >= 0
        accum += I1(s) * c2n(s) * KZCoef(a, q)

    # print(a,c,(-1)**(a+1) / I1(n) * 2 *simplifyWithMaple(mpl,KZCoef(a,0,c)))
    # return (2*simplifyWithMaple(mpl,KZCoef(a,0,c))/I1(n))
    return (-1) ** (a + 1) / I1(n) * (2 * KZCoef(a, 0, c) + accum) * (-1) ** (alpha + beta)  # -1 as alpha+beta is odd

def regZeta(*args):
    args = list(args)
    depth = len(args) -1
    lastArg = args[-1]
    if lastArg == 0:
        if set(args) == {3,5}:
            print("GOT")
        return -zeta(*[x+1 for x in reversed(args[:-1])])
    result = 0
    for i in range(depth):
        newArgs = args.copy()
        newArgs[i] +=1
        newArgs[-1] -= 1
        #print(newArgs)
        result -= Integer(args[i]+1)/Integer(lastArg) * regZeta(*newArgs)
    return result

def modifiedDoubleIntegralPart(n, a, b, c, d):
    result = [Integer(0),Integer(0),Integer(0)]
    for l in range(1, int((a + b + c) / 2) + 1):
        for m in range(1, int((a + b + c) / 2 - l) + 1):
            p = a + b + c - 2 * l - 2 * m
            assert p >= 0
            cof = binom(2 * l, a) * binom(2 * m, a + b - 2 * l) * (Integer(-1)) ** b
            cof -= binom(2 * l, b + c - 2 * m) * binom(2 * m, c) * (Integer(-1)) ** b
            cof += binom(2 * l, b) * binom(2 * m, a) * (Integer(-1)) ** (a + b)
            cof -= binom(2 * l, b) * binom(2 * m, a + b - 2 * l) * (Integer(-1)) ** a
            cof += binom(2 * l, b) * binom(2 * m, c) * (Integer(-1)) ** (c + b)
            cof -= binom(2 * l, b) * binom(2 * m, b + c - 2 * l) * (Integer(-1)) ** c
            if not p==d==0:
                result[p] += cof * J2lm(l, m) * c2n(l) * c2n(m) * KZCoef(p, d)

    for l in range(1, int((b + c + d) / 2) + 1):
        for m in range(1, int((b + c + d) / 2 - l) + 1):
            q = b + c + d - 2 * l - 2 * m
            assert q >= 0 and l > 0 and m > 0
            cof = binom(2 * l, b) * binom(2 * m, b + c - 2 * l) * (Integer(-1)) ** c
            cof += binom(2 * l, c) * binom(2 * m, b) * (Integer(-1)) ** (b + c)
            cof -= binom(2 * l, c) * binom(2 * m, b + c - 2 * l) * (Integer(-1)) ** b
            if not a==q==0:
                result[a] += cof * J2lm(l, m) * c2n(l) * c2n(m) * KZCoef(a, q)

    return result

def rhsDoubleIntegralPart(n,a,b,c,d):
    result = 0
    for l in range(1,int((a+b+c)/2) + 1):
        for m in range(1, int((a + b + c) / 2 - l) + 1):
            p = a+b+c-2*l-2*m
            assert p >= 0
            cof = binom(2*l,a) * binom(2*m,a+b-2*l) * (Integer(-1))**b
            cof -= binom(2*l,b+c-2*m) * binom(2*m,c) * (Integer(-1))**b
            cof += binom(2*l,b) * binom(2*m,a) * (Integer(-1))**(a+b)
            cof -= binom(2*l,b) * binom(2*m,a+b-2*l) * (Integer(-1))**a
            cof += binom(2*l,b) * binom(2*m,c) * (Integer(-1))**(c+b)
            cof -= binom(2*l,b) * binom(2*m,b+c-2*l) * (Integer(-1))**c
            result += cof * J2lm(l,m) * c2n(l) * c2n(m) * KZCoef(p,d)


    for l in range(1, int((b + c + d) / 2) + 1):
        for m in range(1, int((b + c + d) / 2 - l) + 1):
            q = b + c + d - 2 * l - 2 * m
            assert q >= 0
            cof = binom(2*l,b) * binom(2*m,b+c-2*l) * (Integer(-1))**c
            cof += binom(2*l,c) * binom(2*m,b) * (Integer(-1))**(b+c)
            cof -= binom(2*l,c) * binom(2*m,b+c-2*l) * (Integer(-1))**b
            result += cof * J2lm(l, m) * c2n(l) * c2n(m) * KZCoef(a,q)

    return result

def rhsTripleIntegralPart(n,a,b,c,d):
    result = 0
    for h in range(1,(n-3)+1):
        for l in range(1,(n-2-h)+1):
            m = (n-1)-l-h
            assert m >= 1
            cof = binom(2*h,a) * binom(2*l,a+b-2*h) * binom(2*m,d) * (Integer(-1))**(b-d)
            cof += binom(2*h,b) * binom(2*l,b+c-2*h) * binom(2*m,a) * (Integer(-1))**(a+c)
            cof -= binom(2*h,a+d-2*m) * binom(2*l,c) * binom(2*m,d) * (Integer(-1))**(a-c)
            cof += binom(2*h,b) * binom(2*l,a) * binom(2*m,d) * (Integer(-1))**(a+b-d)
            cof -= binom(2*h,b) * binom(2*l,a+b-2*h) * binom(2*m,d) * (Integer(-1))**(a-d)
            cof += binom(2*h,c) * binom(2*l,b) * binom(2*m,a) * (Integer(-1))**(a+b+c)
            cof -= binom(2*h,b) * binom(2*l,a+d-2*m) * binom(2*m,d) * (Integer(-1))**(a+b)
            cof -= binom(2*h,c) * binom(2*l,b+c-2*h) * binom(2*m,a) * (Integer(-1))**(a+b)
            cof += binom(2*h,b) * binom(2*l,c) * binom(2*m,d) * (Integer(-1))**(b+c+d)

            result += cof * K3hlm(h,m,l) * c2n(h) * c2n(l) * c2n(m)

    return result