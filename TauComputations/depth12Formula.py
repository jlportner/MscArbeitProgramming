#from sage.all import *
from zetaLogic import *
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

def KZCoefDepth1(a,b):
    if a+b == 0:
        return 0
    else:
        return (-1)**1 / (2*pi*i)**(a+b+1) * zetaDepth1(a,b)

def KZCoefDepth2(a,b,c):
    return (-1)**2 / (2*pi*i)**(a+b+c+2) * zetaDepth2(a,b,c)

def KZCoefDepth3(a,b,c,d):
    return (-1)**3 / (2*pi*i)**(a+b+c+d+3) * zetaDepth3(a,b,c,d)

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

    return -1 / (I1(n,1) * (2*pi*i)**(2*n+1)) * (binom(2*n,a) + (-1)**(a+1) - binom(2*n,a+1)) * zeta(2*n+1)
