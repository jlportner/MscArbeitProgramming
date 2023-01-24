import numpy as np
import math
from singleValuedKZ import *

def binom(n,k):
    if k < 0 or n < 0:
        return 0
    else:
        return math.comb(n,k)


#n = 10
mpl = startHyperlogProd()


for n in range(4,5):
    #attention a is not 0
    for a in range(1, 2 * n):
        for c in range(0, 2 * n - a):
            b = 2 * n - 1 - c - a

            p1 = 0
            #for i in range(0,n):
            #    p1 += (-1)**(a-c) * (binom(2*n,i+1) - binom(2*n,i) + (-1)**i) * (-1)**(i+1) * (binom(i,a) * binom(2*n-1-i,c) + binom(i,c) * binom(2*n-1-i,a))

            #for i in range(0,b+1):
                #p1 += (-1)**(a-c) * (binom(2*n,a+i+1) - binom(2*n,a+i) + (-1)**(a+i)) * (-1)**(a+i+1) * (binom(a+i,a) * binom(b+c-i,c))
            #    p1 += (-1)**(a-c) * binom(2*n,a+i+1) * (-1)**(a+i+1) * (binom(a+i,a) * binom(b+c-i,c))

            p1 += (-1)**(a+1) * binom(2*n,a) + (-1)**b * binom(2*n,b) + (-1)**(c+1) * binom(2*n,c)

            p2 = 0
            for i in range(0,a+c+1):
                #p2 += binom(a+i,a) * binom(b+c-i,b) * (-1)**(b-i+1) * (binom(2*n,b+c-i) + binom(2*n,a+i) + (-1)**(b+c-i+1)) * (-1/2)
                #p2 += binom(a+i,a) * binom(b+c-i,b) * (-1)**(b-i+1) * (binom(2*n+1,b+c-i+1) + (-1)**(b+c-i+1)) * (-1/2)
                #p2 += (-1)**c / 2 * binom(a+c-i,a) * binom(b+i,i) * (binom(2*n+1,b+i+1) * (-1)**(b+i) - 1)
                p2 += (-1)**(b+c+i) / 2 * binom(a+c-i,a) * binom(b+i,b) * binom(2*n+1,b+i+1)

            p2 -= (-1)**c /2 * binom(2*n,c)
            if p2 != 0:
                print(p1/p2)
            else:
                print(p1,p2)

