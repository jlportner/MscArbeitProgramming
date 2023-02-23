from depth3Formula import *

def doubleIntegralEvenPart(n,a,b,c,d,t=1):
    result = 0
    for s in range(1,n-2 +1):
        l = n-1-s
        assert l >= 1
        for beta in range(s, 2 * s):
            alpha = 2 * s - 1 - beta
            assert 0 <= alpha < beta
            #binom(alpha,) * binom(beta,) * binom(2*l,) * Integer(-1)
            # cab acts on c2l
            interCoef = binom(alpha,a) * binom(beta,a+b-alpha) * binom(2*l,d) * Integer(-1)**(b-d-alpha)
            interCoef -= binom(alpha,a+b-beta) * binom(beta,a) * binom(2*l,d) * Integer(-1)**(b-d-beta)
            interCoef += binom(alpha,b) * binom(beta,b+c-alpha) * binom(2*l,a) * Integer(-1)**(a+c-alpha)
            interCoef -= binom(alpha,a+d-2*l) * binom(beta,c) * binom(2*l,d) * Integer(-1)**(a+c-beta)
            interCoef -= binom(alpha,b+c-beta) * binom(beta,b) * binom(2*l,a) * Integer(-1)**(a+c-beta)
            interCoef += binom(alpha,c) * binom(beta,a+d-2*l) * binom(2*l,d) * Integer(-1)**(a+c-alpha)
            result += J2lm(s, l,t) * cab(s, alpha, beta) * c2n(l) * interCoef * (-1)

            #c2l acts on cab
            interCoef = binom(alpha,a+b-2*l) * binom(beta,d) * binom(2*l,a) * Integer(-1)**(b+d-beta)
            interCoef -= binom(alpha,d) * binom(beta,a+b-2*l) * binom(2*l,a) * Integer(-1)**(b+d-alpha)

            interCoef += binom(alpha,a) * binom(beta,d) * binom(2*l,b) * Integer(-1)**(a+b+d-beta)
            interCoef -= binom(alpha,a+b-2*l) * binom(beta,d) * binom(2*l,b) * Integer(-1)**(a+d-beta)
            interCoef += binom(alpha,a) * binom(beta,a+b-alpha) * binom(2*l,c) * Integer(-1)**(b+c-alpha)
            interCoef -= binom(alpha,a) * binom(beta,d) * binom(2*l,c) * Integer(-1)**(a+c+d-beta)
            interCoef -= binom(alpha,d) * binom(beta,a) * binom(2*l,b) * Integer(-1)**(a+b+d-alpha)
            interCoef += binom(alpha,d) * binom(beta,a+b-2*l) * binom(2*l,b) * Integer(-1)**(a+d-alpha)
            interCoef -= binom(alpha,a+b-beta) * binom(beta,a) * binom(2*l,c) * Integer(-1)**(b+c-beta)
            interCoef += binom(alpha,d) * binom(beta,a) * binom(2*l,c) * Integer(-1)**(a+c+d-alpha)
            result += J2lm(l,s,t) * cab(s,alpha,beta) * c2n(l) * interCoef * (-1) #-1 because alpha+beta is odd

    return result

def ATCoefDepth3(a,b,c,d,mpl=None):
    return interpolationAssociatorDepth3(a,b,c,d,t=Rational(1/2),mpl=mpl)

def interpolationAssociatorDepth3(a,b,c,d,t=1,mpl=None):
    if mpl is None:
        mpl = startHyperlogProd()
    m = a+b+c+d
    if (m+3) % 2 == 1:
        n = int((m+2)/2)
        cabContribution = rhsCabPart(n, a, b, c, d,t)
        cDoubleKZContribution = rhsCDoubleKZpart(n, a, b, c, d,t)
        doubleContribution = rhsDoubleIntegralPart(n, a, b, c, d,t)
        tripleContribution = rhsTripleIntegralPart(n, a, b, c, d,t)
        cabcContribution = cabcPart(n,a,b,c,d,t)
        result = KZCoefDepth3(a,b,c,d) + cabContribution + cDoubleKZContribution + doubleContribution + tripleContribution + cabcContribution
    else:
        n = int((m+3)/2)
        cabContribution = rhsCabPart(n,a,b,c,d,t)
        cDoubleKZContribution = rhsCDoubleKZpart(n,a,b,c,d,t)
        doubleContribution = rhsDoubleIntegralPart(n,a,b,c,d,t)
        doubleEvenContribution = doubleIntegralEvenPart(n,a,b,c,d,t)
        result = (KZCoefDepth3(a,b,c,d)  + cabContribution + cDoubleKZContribution + doubleContribution + doubleEvenContribution)
    result = simplifyWithMaple(mpl, result)
    result = convertToBrownBasis(result)
    return result

if __name__ == "__main__":
    print(ATCoefDepth3(0,0,1,4))
    exit()
    mpl = startHyperlogProd()
    for n in range(6,7):
        m = n - 3
        print("\\subsection{words of length %d}" % n)
        print("\\begin{align*}")
        for a in range(0, m + 1):
            for b in range(0, m + 1 - a):
                for d in range(0, m + 1 - a - b):
                    c = m - a - b - d
                    x = interpolationAssociatorDepth3(a, b, c, d,t=Rational(1/2),mpl=mpl)
                    latexStr = latex(x)
                    word = ""
                    for i in (a,b,c,d):
                        if i == 0:
                            word += "y "
                        else:
                            word += "x^%d y " % i
                    word = word[:-2].replace("y y y","y^3").replace("y y","y^2")
                    print("\Phi_{AT} \mid_{%s}" % word, "&=", latexStr, "\\\\")
        print("\\end{align*}")
        
    exit()
    #"""

    mpl = startHyperlogProd()
    n = 9
    m = n-3
    #a,b,c,d = 4,2,0,4
    for a in range(0,m+1):
        for b in range(0,m+1-a):
            for d in range(0,m+1-a-b):
                c = m - a - b - d
                ogD = rhsDoubleIntegralPart(n,a,b,c,d,t=1)
                atT = rhsDoubleIntegralPart(n,a,b,c,d,t=Rational(1/2))
                atD = rhsTripleIntegralPart(n,a,b,c,d,t=1)
                ogT = rhsTripleIntegralPart(n,a,b,c,d,t=Rational(1/2))

                if ogD != 0 and ogT != 0:
                    print(atD/ogD, atT/ogT)
                else:
                    if ogD == 0 and ogT != 0:
                        print("Double og zero, at :", atD, "Triple:", atT/ogT)
                    if atT == 0 and atD != 0:
                        print("Double:", atD/ogD, "Triple og zero, at:", atT)
                    else:
                        print("all zero")

                #C1 = interpolationAssociatorDepth3(a,b,c,d,mpl=mpl)
                #C2 = convertToBrownBasis(simplifyWithMaple(mpl,Integer(-1)**(a+b+c+d+3) * KZCoefDepth3(a,b,c,d)))
                #print(C1)
                #print(C2)
                #exit()
                #print(C1 - C2)
    exit()