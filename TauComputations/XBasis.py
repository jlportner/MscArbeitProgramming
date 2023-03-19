from depth3Formula import *

class LieElement(dict):
    #expr = {}
    def __init__(self,data=None):
        if data is not None:
            super().__init__(data)
        else:
            super().__init__()

    @staticmethod
    def x():
        return LieElement({"x":1})

    @staticmethod
    def y():
        return LieElement({"y":1})

    def __copy__(self):
        #cls = self.__class__
        #result = cls.__new__(cls)
        #result.__dict__.update(self.__dict__)
        #super.copy()
        return LieElement(self)

    def copy(self):
        return self.__copy__()

    def __iadd__(self, other):
        if isinstance(other,self.__class__):
            for key in other:
                if key in self:
                    self[key] += other[key]
                else:
                    self[key] = other[key]
            return self
        else:
            return NotImplemented

    def __add__(self, other):
        A = self.copy()
        A += other
        return A

    def __isub__(self, other):
        if isinstance(other,self.__class__):
            for key in other:
                if key in self:
                    self[key] -= other[key]
                else:
                    self[key] = -other[key]
            return self
        else:
            return NotImplemented

    def __sub__(self, other):
        A = self.copy()
        A -= other
        return A

    def __neg__(self):
        return -1 * self


    def __imul__(self, other):
        if isinstance(other, self.__class__):
            return NotImplemented
        for x in self:
            self[x] *= other
        return self

    def __mul__(self, other):
        if isinstance(other,self.__class__):
            result = {}
            for y in other:
                for x in self:
                    result += LieElement({x+y:other[y]*self[x]})
        else:
            result = self.copy()
            result *= other

        return result

    def __rmul__(self, other):
        return self * other

def lieAd(s):
    if s == 0:
        return LieElement.y()
    else:
        expr = lieBracket(LieElement.x(),LieElement.y())
        for i in range(s-1):
            expr = lieBracket(LieElement.x(),expr)

    return expr

def lieBracket(A,B):
    result = LieElement()
    for x in A:
        for y in B:
            result += LieElement({x+y: A[x] * B[y]})
            result += LieElement({y+x:-A[x] * B[y]})
    return result

def assAction(A,B):
    result = LieElement()
    for u in A:
        for v in B:
            result += LieElement({u+v:A[u]*B[v]})
            for i in range(len(v)):
                if v[i] == "y":
                    result += LieElement({v[:i+1] + u + v[i+1:] : A[u] * B[v],v[:i] + u + v[i:] : -A[u] * B[v]})
    return result

def IharaBracket(A,B):
    return assAction(A,B) - assAction(B,A)

def adjustCabc(alpha,beta,gamma,n):
    res = cabc(alpha,beta,gamma)
    res = convertToSingleValuedSchnetz(res)
    if n == 5:
        e1,e2 = var("e1,e2")
        res = res.subs(zetaSV(5,3,3) == e2 - Rational(22020/3553) * zeta(3)**2 * zeta(5), zeta(11)==e1)
        res = res.expand()
        res = res.subs(e1==0)
    elif n==6:
        e1,e2,e3 = var("e1,e2,e3")
        res = res.subs(zeta(13)==e1,zetaSV(5,5,3) == e2 + Rational(203950/5681) * zeta(5)**2 * zeta(3))
        res = res.subs(zetaSV(7,3,3) == e3 + Rational(244740/5681) * zeta(5)**2 * zeta(3) - Rational(123508/7429)* zeta(7) * zeta(3)**2)
        res = res.expand()
        res = res.subs(e1 == 0)
    return res

def bracketingStrToElement(expr):
    opener = -1
    comma = -1
    depth = 0
    for i in range(len(expr)):
        if expr[i] == '[':
            if opener == -1:
                opener = i
            else:
                depth +=1
        if expr[i] == ']':
            if depth == 0:
                if comma == -1:
                    raise ValueError("Missing ',' in expression")
                else:
                    A = bracketingStrToElement(expr[opener+1:comma])
                    B = bracketingStrToElement(expr[comma+1:i])
                    return lieBracket(A,B)
            else:
                depth -= 1

        if opener != -1 and expr[i] == ',' and depth == 0:
            comma = i

    return LieElement({expr.strip():1})

#END of Lie Element functions

def tauDepth3(n):
    m = 2*n-2
    result = LieElement()
    for gamma in range(int(m / 3 + 1), m + 1):
        for beta in range(min(m - gamma + 1, gamma)-1,max(0,m-2*gamma)-1,-1):
            alpha = m - beta - gamma
            assert 0 <= beta < gamma and 0 <= alpha <= gamma
            expr = lieBracket(lieAd(alpha),lieBracket(lieAd(beta),lieAd(gamma)))
            result += adjustCabc(alpha,beta,gamma,n) *expr

    return result

def adjustCabcForSigma(alpha,beta,gamma,n):
    res = cabc(alpha,beta,gamma)
    res = convertToSingleValuedSchnetz(res)
    if n==4:
        e1= var("e1")
        res = res.subs(zeta(9)==e1 * (2*pi*i)**9)
        res = res.expand()
    elif n == 5:
        e1,e2 = var("e1,e2")
        res = res.subs(zetaSV(5,3,3) == e2 - Rational(22020/3553) * zeta(3)**2 * zeta(5), zeta(11)==e1 * (2*pi*i)**11)
        res = res.expand()
        res = res.subs(e2==0)
    elif n==6:
        e1,e2,e3 = var("e1,e2,e3")
        res = res.subs(zeta(13)==e1 * (2*pi*i)**13,zetaSV(5,5,3) == e2 + Rational(203950/5681) * zeta(5)**2 * zeta(3))
        res = res.subs(zetaSV(7,3,3) == e3 + Rational(244740/5681) * zeta(5)**2 * zeta(3) - Rational(123508/7429)* zeta(7) * zeta(3)**2)
        res = res.expand()
        res = res.subs(e2==0,e3==0)
    return res
def adjustCabForSigma(n,alpha,beta):
    res = cab(n,alpha,beta)
    if n == 5:
        e1 = var("e1")
        res = res.subs(zeta(11)==e1 * (2*pi*i)**11)
        res = res.expand()
    elif n==6:
        e1 = var("e1")
        res = res.subs(zeta(13)==e1 * (2*pi*i)**13)
        res = res.expand()
    return res

def sigmaNContrib(n):

    normalizer = c2n(n)
    e1 = var("e1")
    normalizer = normalizer.subs(zeta(2*n+1) == e1 * (2*pi*i)**(2*n+1))
    print(normalizer)
    print(latex(normalizer/ normalizer), "\\ad_x^{%d}(y)" % (2*n))

    m = 2*n-1
    for beta in range(m, n-1,-1):
        alpha = m-beta
        assert 0 <= alpha < beta
        expr = "[\\ad_x^{%d}(y),\\ad_x^{%d}(y)]" % (alpha, beta)
        print(latex(adjustCabForSigma(n,alpha,beta)/normalizer),expr)

    m = 2 * n - 2
    for gamma in reversed(range(int(m / 3 + 1), m + 1)):
        for beta in range(min(m - gamma + 1, gamma) - 1, max(0, m - 2 * gamma) - 1, -1):
            alpha = m - beta - gamma
            assert 0 <= beta < gamma and 0 <= alpha <= gamma
            expr = "[\\ad_x^{%d}(y),[\\ad_x^{%d}(y),\\ad_x^{%d}(y)]]" % (alpha,beta,gamma)
            print(latex(adjustCabcForSigma(alpha,beta,gamma,n) / normalizer),expr)

def sigma5Depth3():
    return 10*bracketingStrToElement("[x,[x,[x,[x,y]]]]") + 5 * bracketingStrToElement("[[x,y],[x,[x,y]]]") + 20 * bracketingStrToElement("[y,[x,[x,[x,y]]]]")

def sigma3Depth3():
    return -24*bracketingStrToElement("[x,[x,y]]") + 24 * bracketingStrToElement("[[x,y],y]")

def sigma7Depth3():
    return -14 * bracketingStrToElement("[x, [x, [x, [x, [x, [x, y ]]]]]]") - 42 * bracketingStrToElement("[[x,[x,y]],[x,[x,[x,y]]]]") \
        - 56 * bracketingStrToElement("[[x,y],[x,[x,[x,[x,y]]]]]") - 42 * bracketingStrToElement("[y,[x,[x,[x,[x,[x,y]]]]]]")

def discardHigherDepths(A,n):
    result = {}
    for x in A:
        if x.count("y") <= n:
            result[x] = A[x]

    return LieElement(result)

def discardZeros(A):
    return LieElement({key:A[key] for key in A if A[key] != 0})

def checkCorrectGeneratorsFor11():
    e1, e2 = var("e1,e2")
    s3 = sigma3Depth3()
    s5 = sigma5Depth3()
    I1 = IharaBracket(s3, s5)
    I1 = discardHigherDepths(I1, 2)
    s335 = IharaBracket(s3, I1)
    s335 = discardHigherDepths(s335, 3)
    s335 = discardZeros(s335)
    s335 *= Rational(-323323/2400) * e2 / (2*pi*i)**11
    X5 = tauDepth3(5)
    X5 = discardZeros(X5)

    temp = X5 - s335
    temp = discardZeros(temp)
    print(temp)

def checkCorrectGeneratorsFor13():
    e1, e2, e3 = var("e1,e2,e3")

    s3 = sigma3Depth3()
    s5 = sigma5Depth3()
    s7 = sigma7Depth3()
    I1 = IharaBracket(s5, s3)
    I1 = discardHigherDepths(I1, 2)
    s553 = IharaBracket(s5, I1)
    s553 = discardHigherDepths(s553, 3)
    s553 = discardZeros(s553)
    s553 *= Rational(-676039/600)*e2/(2*i*pi)**13

    I1 = IharaBracket(s3, s7)
    I1 = discardHigherDepths(I1, 2)
    s733 = IharaBracket(s3, I1)
    s733 = discardHigherDepths(s733, 3)
    s733 = discardZeros(s733)
    s733 *= -Rational(482885/672)*e2/(2*i*pi)**13 + Rational(2414425/4032)*e3/(2*i*pi)**13

    X13 = tauDepth3(6)
    X13 = discardZeros(X13)
    temp = X13 - s553 - s733
    temp = discardZeros(temp)
    print(temp)

