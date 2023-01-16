import numpy as np
import math

def binom(n,k):
    if k < 0 or n < 0:
        return 0
    else:
        return math.comb(n,k)


#n = 10



for n in range(1,10):
    for a in range(n):
        c = 2 * n - 1 - a
        assert(c > a)
        #r1 = binom(2*s,c) + (-1)**a * binom(2*s,a)**2
        #r2 = binom(2*s,a) * binom(2*n-2 -2*s,c) * (-1)**(a+1)
        #print(r1,r2,r2-r1)
        result = 0
        for k in range(c + 1):
            result += binom(2 * n - (k + 1), a) * binom(2 * n + 1, k + 1) * (-1) ** k

        print(result - binom(2*n,a) + (-1)**a)
        
print(result)