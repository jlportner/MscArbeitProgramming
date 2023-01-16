import numpy as np
from sage.all import *
import sage.modular.multiple_zeta as mzv
from mapleHandler import *
from timeit import default_timer as timer
from zetaLogic import *
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

#w = word.Word("xyx.wo")
#mzv.Multizeta(w)



def regZetaSchnetz(*args):
    args = list(args)
    depth = len(args) -1
    lastArg = args[-1]
    if lastArg == 0:
        return zeta(*[x+1 for x in args[:-1]])
    result = 0
    for i in range(depth):
        newArgs = args.copy()
        newArgs[i] +=1
        newArgs[-1] -= 1
        #print(newArgs)
        result -= Integer(args[i]+1)/Integer(lastArg) * regZetaSchnetz(*newArgs)
    return result

def regZetaWithHyperlog(*args,mpl=None):
    args = list(args)
    length = sum(args)
    depth = len(args) -1
    xyWord = np.zeros(length+depth+2,dtype=int)
    xyWord[0] = 1
    for i in range(depth):
        index = 1+ sum(args[:i+1]) + i
        xyWord[index] = 1

    #print(xyWord)
    #xyWord = np.flip(xyWord)
    expr = mpl("Convert(i" + np.array2string(xyWord,separator=",") + ",zeta)")
    return (-1)**depth * parseMapleToSage(str(expr))
    #return mzv.All_iterated(QQ)(tuple(xyWord)).regularise().composition()


def mzeta(*args):
    x = [Integer(i) for i in args]
    return mzv.Multizetas(SR)(tuple(x))

def KZCoefDepth1(a,b):
    if a+b == 0:
        return 0
    else:
        return (-1)**1 / (2*pi*i)**(a+b+1) * zetaDepth1(a,b)

def KZCoefDepth2(a,b,c):
    return (-1)**2 / (2*pi*i)**(a+b+c+2) * zetaDepth2(a,b,c)

def KZCoefDepth3(a,b,c,d):
    return (-1)**3 / (2*pi*i)**(a+b+c+d+3) * zetaDepth3(a,b,c,d)

"""
def KZCoef(*args):
    args = list(args)
    numberOfx = sum(args)
    depth = len(args) - 1
    if depth + numberOfx < 2:
        return 0
    else:
        return (-1)**depth / (2 * pi * i)**(numberOfx + depth) * regZetaSchnetz(*args)
"""
def I1(n,t=1):
    if t==1:
        return Integer(factorial(2*n)**2) / Integer(factorial(4*n+1))
    if t == Rational(1/2):
        return Integer(factorial(2*n)**2) / (Integer(2)*Integer(factorial(4*n+1)))
    else:
        x = var("x")
        assume(x >= 0, x <= t)
        f = (x * (x - 1)) ** (2 * n)
        return integral(f, x, 0, t)

def J2lm(l,m,t=1):
    x,y = var("x,y")
    assume(x>=0,x<=t)
    f = (x*(x-1))**(2*l) * (y*(y-1))**(2*m)
    return integral(integral(f,y,0,x),x,0,t)

def K3hlm(l,m,h,t=1):
    x,y,z = var("x,y,z")
    assume(x>=0,x<=t)
    assume(y>=0,y<=t)
    f = (x*(x-1))**(2*l) * (y*(y-1))**(2*m) * (z*(z-1))**(2*h)
    return integral(integral(integral(f,z,0,y), y, 0, x), x, 0, t)

def c2n(n):
    assert n > 0
    return 2*(4*n+1) * binom(4*n,2*n) / (2*pi*i)**(2*n+1) * zeta(2*n+1)

def cab(n,alpha,beta = -1):
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
    for k in range(c+1):
        accum += binom(2*n-(k+1),a) * (binom(2*n+1,k+1) * (-1)**k -1)

    return -1 / (I1(n,1) * (2*pi*i)**(2*n+1)) * accum * zeta(2*n+1)

def cabGuess(n,alpha,beta = -1):
    assert n > 0
    if beta == -1:
        beta = 2 * n - 1 - alpha
    elif beta + alpha != 2 * n - 1:
        raise ValueError("alpha + beta = 2n-1 needs to be satisfied!")
    if alpha >= beta:
        raise ValueError("alpha < beta needs to be satisfied!")
    a = alpha

    return -1 / (I1(n,1) * (2*pi*i)**(2*n+1)) * (binom(2*n,a) + (-1)**(a+1) - binom(2*n,a+1)) * zeta(2*n+1)

if __name__ == "__main__":
    #s = timer()
    #print(regZeta(1,1,18))
    #print(mzeta(5,3,3,4,2))
    #e = timer()
    #print("Took: " + str(e - s) + "s")
    #exit()

    #mpl = startHyperlogProd()
    #tpl = (0,4,3,3)
    #print(simplifyWithMaple(mpl,zetaDepth3(*tpl)))
    #print(simplifyWithMaple(mpl,regZetaSchnetz(*tpl)))
    #print(simplifyWithMaple(mpl,regZetaWithHyperlog(*tpl)))


    for n in range(1,20):
        print("\\subsection{n=%d}" % n,)
        print("\\begin{align*}")
        for beta in range(n, 2 * n):
            alpha = 2 * n - 1 - beta
            ck = cab(n,alpha, beta)
            cGuess= cabGuess(n, alpha, beta)
            print(ck-cGuess)
            #print("c_{%d,%d} &=" % (alpha,beta),latex(ck), "\\\\")
        print("\\end{align*}")

