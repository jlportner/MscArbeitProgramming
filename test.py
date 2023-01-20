import numpy as np
import math

def binom(n,k):
    if k < 0 or n < 0:
        return 0
    else:
        return math.comb(n,k)


#n = 10



for n in range(4,5):
    #attention a is not 0
    for a in range(1, 2 * n):
        for c in range(0, 2 * n - a):
            b = 2 * n - 1 - c - a

            p1 = 0
            for i in range(0,n):
                p1 += (-1)**(a-c) * (binom(2*n,i+1) - binom(2*n,i) + (-1)**i) * (-1)**(i+1) * (binom(i,a) * binom(2*n-1-i,c) + binom(i,c) * binom(2*n-1-i,a))

            p2 = 0
            for i in range(0,c+1):
                p2 += binom(a+i,a) * binom(b+c-i,b) * (binom(2*n,b+c-i) + binom(2*n,a+i) + (-1)**a)

            if p2 != 0:
                print(p1/p2)
            else:
                print(p1,p2)