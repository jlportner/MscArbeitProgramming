import numpy as np
from depth3Formula import *
from itertools import combinations

def rhsDoubleTripleDiff(n,a,b,c,d):
    result = 0
    if (a,b,c,d) == (0,5,2,1):
        print("")
    for l in range(1,int((a+b+c)/2.0)):
        for m in range(1,int((a+b+c)/2.0) - l + 1):
            p = a+b+c-2*l-2*m
            assert p >= 0
            interCoef = binom(2*l,c) * binom(2*m,b+c-2*l) * Integer(-1)**b
            if p + d != 0:
                result += (-1)**d * binom(p+d,p)  * interCoef

    for l in range(1,(n-3)+1):
        for m in range(1,(n-2-l)+1):
            h = (n-1)-m-l
            assert h >= 1
            cof = binom(2*l,c) * binom(2*m,a+d-2*h) * binom(2*h,d) * Integer(-1)**(a-c)
            result -= cof
    return result

def checkFullMat(n,mpl):
    M,rhs = fullMatAndRhs(n,mpl)
    #x = cDepth3(n,mpl)
    A = matrix(SR,M)
    #print(A)
    b = matrix(SR,rhs).T
    #Ab = A.augment(b,subdivide=true)
    #return Ab#Ab.echelon_form()
    return A.solve_right(b)

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start

def weight4Check(mpl=None):
    n=4
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
                mL += [rhsTripleIntegralPartNew(n,a,b,c,d)]
                mL += modifiedDoubleIntegralPart(n,a,b,c,d)
                line = rhsLine(n,a,b,c,d)
                line = simplifyWithMaple(mpl, line)
                M += [mL]
                rhs += [line]
                #print(mL, line)

    #M = np.array(M, dtype=int)

    A = matrix(M)
    print(A)
    b = matrix(rhs).T
    print(b)
    return A.solve_right(b)


#mpl = startHyperlogProd()
#M, rhs = fullMatAndRhs(4, mpl)
#data = np.hstack((np.array(M),np.array(rhs).reshape(84,1)))
#np.save("n4Data",data,allow_pickle=True)
#exit()

def findNotWorkingLines(mpl=None):
    data = np.load("n4Data.npy",allow_pickle=True)
    M = data[:,:-1].tolist()
    rhs = data[:,-1].tolist()

    A = matrix(SR, M)
    b = matrix(SR, rhs).T
    for ind in combinations(range(84),3):
        ASub = A.matrix_from_rows(ind)
        bSub = b.matrix_from_rows(ind)
        if ASub.rank() != 3:
            try:
                A.solve_right(b)
            except:
                print(ASub,bSub)
                #return

"""
for a in range(0,m+1):
    for b in range(0,m+1-a):
        for d in range(0,m+1-a-b):
            c = m - a - b - d
            assert c >= 0 and a+b+c+d == m
            #if (a,b,c,d) == (0,6,0,0):
            #    breakpoint()
            tri = rhsTripleIntegralPartNew(n,a,b,c,d)
            tri = simplifyWithMaple(mpl,tri)
            dou = rhsDoubleIntegralPartNewest(n,a,b,c,d)
            dou = simplifyWithMaple(mpl,dou)
            dif = (tri + dou) + Rational(1/2) * tri
            #line = rhsLine(n,a,b,c,d)
            #line = simplifyWithMaple(mpl,line)
            #print(a,b,c,d,line,tri)
            print(a,b,c,d,dif)

exit()
"""

"""
for a in range(0,m+1):
    for b in range(0,m+1-a):
        for d in range(0,m+1-a-b):
            c = m - a - b - d
            assert c >= 0 and a+b+c+d == m
            D1 = rhsDoubleTripleDiff(n,a,b,c,d)
            print(a,b,c,d,D1)

exit()
"""

#print(weight4Check())
#print(checkFullMat(n,mpl))
#exit()

#good testing candidates
# 0 6 0 0 and 0 0 6 0
# 5 1 0 0 and 0 0 1 5
"""

f1 = simplifyWithMaple(mpl,rhsLine(n,0,6,0,0))
f2 = simplifyWithMaple(mpl,rhsLine(n,0,0,6,0))
print(f1,"\n",f2)



m1=np.array(matLine(m,0,0,6,0))
m2=np.array(matLine(m,0,6,0,0))
print(m1,m2)

print(m1 - m2)
print(f1 - f2)
exit()

"""
#print(rhsLine(n,4,0,0,4))
#exit()
mpl = startHyperlogProd()
def adjustToBasis(alpha,beta,gamma, n):
        res = cabc(alpha, beta, gamma)
        res = convertToSingleValuedSchnetz(res)
        if n == 4:
            e1 = var("e1")
            res = res.subs(zeta(9) == e1 * (2 * pi * i) ** 9)
            res = res.expand()
        elif n == 5:
            e1, e2 = var("e1,e2")
            res = res * (2 *pi *I)**11
            res = res.subs(zetaSV(5, 3, 3) == (Rational(5/7759752) * e2 - Rational(22020 / 3553) * zeta(3) ** 2 * zeta(5)),
                           zeta(11) == Rational(1/116396280) * e1)
            res = res.expand()
            #res = res.subs(e2 == 0)
        elif n == 6:
            e1, e2, e3 = var("f1,f2,f3")
            res = res * (2 * pi * I) ** 13
            res = res.subs(zeta(13) == Rational(1/2974571600) * e1,
                           zetaSV(5, 5, 3) == Rational(1/13520780) * e2 + Rational(203950 / 5681) * zeta(5) ** 2 * zeta(3))
            res = res.subs(zetaSV(7, 3, 3) == Rational(1/19315400) * e3 + Rational(244740 / 5681) * zeta(5) ** 2 * zeta(3) - Rational(
                123508 / 7429) * zeta(7) * zeta(3) ** 2)
            res = res.expand()
            #res = res.subs(e2 == 0, e3 == 0)
        elif n == 7:
            e1, e2, e3 = var("g1,g2,g3")
            res = res * (2 * pi * I) ** 15
            res = res.subs(zeta(15) == Rational(1/70578471600) * e1,
                           zetaSV(9,3,3) == Rational(1/258529200) * e2 + 48 * zeta(5)**3 + Rational(8694314/37145) * zeta(7) * zeta(5) * zeta(3) - Rational(68094/2185) * zeta(9) * zeta(3)**2)
            res = res.subs(zetaSV(7,3,5) == Rational(7/1706292720) * e3 - 56 * zeta(5)**3 - Rational(5826772/22287) * zeta(7) * zeta(5) * zeta(3))
            res = res.expand()
        return res

#createCache()
for n in range(5,6):
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

"""

for n in range(5, 6):
    #n=2
    print("\\subsection{n=%d}" % n, )
    print("\\begin{align*}")
    sol,varList = cDepth3(n,mpl=mpl)
    for j,x in enumerate(sol):
        x = convertToSingleValuedSchnetz(x)
        x = expand(x * (2 *pi * I)**(2*n+1))
        latexStr = latex(x)
        if n >= 5:
            i=1
            termsPerLine = 3
            if n >= 7:
                termsPerLine = 3
            breakInd = find_nth(latexStr, "\\frac", termsPerLine*i+1) - 2
            while breakInd > 0:
                latexStr = latexStr[:breakInd] + "\\\\ &" + latexStr[breakInd:]
                i+=1
                breakInd = find_nth(latexStr, "\\frac", termsPerLine * i + 1) - 2
        print("c_{%d,%d,%d}" % varList[j], "&=", latexStr, "\\\\")
    print("\\end{align*}")
"""
#print(sol)
exit()
n=4
m = 2 * n - 2
M = fullRankMat(m)
print(M)
print(np.linalg.matrix_rank(M,tol=0),M.shape)
exit()