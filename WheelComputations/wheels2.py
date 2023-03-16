import math

def binom(n,k):
    if k < 0:
        return 0
    else:
        return math.comb(n,k)



def addToDict(dic,w,coef):
    if w in dic:
        dic[w] += coef
    else:
        dic[w] = coef
    #return dic

def propF(C):
    C = list(C.items())
    result = {}
    for el, coef in C:
        a = el[0]
        z = el[1]
        word = el[2:]
        if z == "z":
            e1 = "0c" + word + a
            e2 = a + "z0" + word
            addToDict(result,e1,-coef)
            addToDict(result,e2,coef)
        if z == "c":
            e1 = "0z" + a + word
            e2 = a + "c" + word + "0"
            addToDict(result,e1,coef)
            addToDict(result,e2,-coef)
    return result

def f(n):
    C = {"1z": 1, "1c": -1}
    for i in range(n):
        C = propF(C)
    return C

def claimf(k):
    C = {}
    for n in range(k):
        m = k-1-n
        addToDict(C,"0z" + "0" * n + "1" + "0" * m, (-1)**(m+1) * binom(k,n))
        addToDict(C,"0c" + "0" * n + "1" + "0" * m, (-1)**(m+1) * binom(k, m))
    return C | {"1z" + "0" * k: 1, "1c" + "0" * k: (-1)**(k+1)}

def getNMFromWord(w):
    ind = w.find("1")
    return ind, len(w) - ind - 1

def propFTil(C):
    C = list(C.items())
    result = {}
    for el, coef in C:
        if len(el) > 1:
            a = el[0]
            z = el[1]
        else:
            z=""

        if z == "z":
            word = el[2:]
            #edge to next
            addToDict(result, "0" + word + a, -coef)
            #edge to center
            if a == "0":
                n,m = getNMFromWord(word)
                addToDict(result, "0" + word + "0", (n+1) * -coef)
                addToDict(result, word + "00", (m+1) * -coef)
            else:
                n = len(word)
                addToDict(result, word + "01", (n+1) * -coef)
        elif z == "c":
            word = el[2:]
            #edge to next
            addToDict(result, a + word + "0", -coef)
            #edge to center
            if a == "0":
                n,m = getNMFromWord(word)
                addToDict(result, "00" + word, (n+1) * -coef)
                addToDict(result, "0" + word + "0", (m+1) * -coef)
            else:
                n = len(word)
                addToDict(result, "10" + word, (n+1) * -coef)
        else:
            word = el
            addToDict(result, word + "0", -coef) #100% correct
            addToDict(result, "0" + word, coef) #100% correct
    return result

def fTilO(n):
    C = {"1": -1}
    for i in range(n-1):
        fn = f(i)
        for e in fn:
            addToDict(C,e, (-1)**i * fn[e])
        C = propFTil(C)
    return C

def fTil(n):
    C = {"1": -1}
    for i in range(n-1):
        fn = f(i)
        for e in fn:
            addToDict(C,e, fn[e])
        C = propFTil(C)
        #for e in C:
        #    C[e] *= -11
    return C


def claimfTil(k):
    result = {}
    for n in range(0,k):
        m = k-1-n
        result[n*"0" + "1" + m * "0"] = binom(n+m,n) * (-1)**(m+1) * (2*k-1)
    return result

def wheelValue(k):
    C = fTil(k)
    result = 0
    for e in C:
        n,m = getNMFromWord(e)
        result += C[e] * (-1)**m * binom(n+m,n)
    return result
if __name__ == "__main__":
    #print(f(1))
    #print(fTil(2))
    print(f(5))
    print(claimf(5))
    exit()
    n=3
    print(fTil(n))
    print(claimfTil(n))
    print(wheelValue(n))
    print("")
    #print(res)
    #exit()