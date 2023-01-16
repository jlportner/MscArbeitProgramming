from sage.all import *
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
"""

import sage.modular.multiple_zeta as mzv

bDeg1 = []
bDeg2 = [mzv.Multizeta(2)]
bDeg3 = [mzv.Multizeta(3)]
bDeg4 = [mzv.Multizeta(2) ** 2]
bDeg5 = [mzv.Multizeta(5), mzv.Multizeta(3) * mzv.Multizeta(2)]
bDeg6 = [mzv.Multizeta(3) ** 2, mzv.Multizeta(2) ** 3]
bDeg7 = [mzv.Multizeta(7), mzv.Multizeta(5) * mzv.Multizeta(2), mzv.Multizeta(3) * mzv.Multizeta(2) ** 2]
bDeg8 = [mzv.Multizeta(3,5),mzv.Multizeta(3,5),mzv.Multizeta(3)**2 * mzv.Multizeta(2), mzv.Multizeta(2)**2]
bDeg9 = []
bDeg10= []
base = [bDeg1, bDeg2, bDeg3, bDeg4, bDeg5, bDeg6, bDeg7, bDeg8, bDeg9, bDeg10]

def brownBasis(n):
    assert n <= 10
    return base[n-1]
"""
def parseMapleToSage(inStr):
    inStr = list(inStr)
    i = 3
    while i < len(inStr):
        if inStr[i-3:i+1] == list("zeta"):
            if inStr[i +1] != "[":
                raise ValueError("Error: Not proper expression")
            else:
                inStr[i+1] = "("
            brackCount = 1
            j = i+1
            while brackCount != 0 and j < len(inStr):
                j+=1
                if inStr[j] == '[':
                    brackCount +=1
                elif inStr[j] == ']':
                    brackCount -= 1
            if brackCount == 0:
                inStr[j] = ')'
            else:
                raise ValueError("Error: Not proper expression 2")
            i+=4
        if inStr[i-1:i+1] == list("Pi"):
            inStr[i-1] = 'p'
        i+=1
    return SR(''.join(inStr))

def parseStringToMaple(inStr):
    inStr = list(inStr)
    i = 3
    while i < len(inStr):
        if inStr[i-3:i+1] == list("zeta"):
            if inStr[i +1] != "(":
                raise ValueError("Error: Not proper expression")
            else:
                inStr[i+1] = "["
            brackCount = 1
            j = i+1
            while brackCount != 0 and j < len(inStr):
                j+=1
                if inStr[j] == '(':
                    brackCount +=1
                elif inStr[j] == ')':
                    brackCount -= 1
            if brackCount == 0:
                inStr[j] = ']'
            else:
                raise ValueError("Error: Not proper expression 2")
            i+=4
        if inStr[i-1:i+1] == list("pi"):
            inStr[i-1] = 'P'
        i+=1
    return ''.join(inStr)



def simplifyWithMaple(mpl,expr):
    exStr = str(expr)
    exStr = parseStringToMaple(exStr)

    optimizedExpr = mpl("Convert(Convert(" + exStr + ",f),zeta)")
    return parseMapleToSage(str(optimizedExpr))

def startHyperlogProd():
    mpl = Maple(script_subdirectory="/media/Transfer/ETH/MscProgramming/hyperlog/")
    mpl.read("HyperlogProcedures")
    return mpl
