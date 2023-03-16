import networkx as nx
import numpy as np
import sympy as sp
from GCboundary import *
from tderIdentification import *
path = "/media/Transfer/ETH/MscArbeit/Images/"

def readGraphfromEdgeList(edges):
    eL = []
    for i,e in enumerate(edges):
        eL += [str(e[0]) + " " + str(e[1]) + " {'order': " + str(i+1) + "}"]

    return nx.parse_edgelist(eL,create_using=nx.MultiGraph,nodetype=int)




inp = "16 17 18 23 25 28 34 38 46 48 57 58 68 78 1 12 13 18 25 26 37 38 45 46 47 56 57 68 78 -7 12 14 18 23 27 35 37 46 48 57 58 67 68 78 -21/8 12 14 16 23 25 36 37 45 48 57 58 67 68 78 77/8 13 14 18 23 25 28 37 46 48 56 57 67 68 78 -77/4 13 16 17 24 25 26 35 37 45 48 58 67 68 78 -7 12 13 15 24 27 35 36 46 48 57 58 67 68 78 -35/8 14 15 17 23 26 28 37 38 46 48 56 57 68 78 49/4 12 13 18 24 26 37 38 46 47 56 57 58 68 78 49/8 12 16 18 27 28 34 36 38 46 47 56 57 58 78 -147/8 14 17 18 23 25 26 35 37 46 48 56 58 67 78 77/8 12 15 16 27 28 35 36 38 45 46 47 57 68 78 -21/8 12 13 18 26 27 35 38 45 46 47 56 57 68 78 -105/8 12 14 18 23 27 35 36 45 46 57 58 67 68 78 -35/8 12 14 18 23 27 36 38 46 48 56 57 58 67 78 7/8 14 15 16 23 26 28 37 38 46 48 57 58 67 78 -49/4 12 14 15 23 27 35 36 46 48 57 58 67 68 78 35/8 12 15 18 23 28 34 37 46 48 56 57 67 68 78 105/8 12 13 14 27 28 36 38 46 47 56 57 58 68 78 -49/8 12 14 17 23 26 37 38 46 48 56 57 58 68 78 -49/8 12 13 18 25 27 34 36 47 48 56 58 67 68 78 35/4 12 16 18 25 27 35 36 37 45 46 48 57 68 78 49/16 12 13 14 25 26 36 38 45 47 57 58 67 68 78 -119/16 12 13 18 25 27 35 36 46 47 48 56 57 68 78 7 12 13 15 24 28 36 38 47 48 56 57 67 68 78 49/8 12 14 18 25 28 34 36 38 47 57 58 67 68 78 -7 12 13 14 23 28 37 46 48 56 57 58 67 68 78 77/4 12 16 18 25 27 35 36 37 45 46 48 58 67 78 -77/16 12 15 17 25 26 35 36 38 45 47 48 67 68 78 -49/8 12 14 18 23 27 35 38 46 47 57 58 67 68 78 77/4 13 15 18 24 26 28 37 38 46 47 56 57 68 78 -49/4 12 14 15 23 27 36 38 46 48 57 58 67 68 78 35/2 13 14 18 25 26 28 36 38 47 48 56 57 67 78 -49/4 12 13 18 25 27 34 36 46 48 57 58 67 68 78 -105/8 12 14 18 23 28 35 37 46 48 56 57 67 68 78 -7 12 15 16 25 27 35 36 38 46 47 48 57 68 78 -7 12 14 18 23 28 36 38 46 47 56 57 58 67 78 -7 12 13 16 25 28 34 37 47 48 57 58 67 68 78 -147/16 12 15 16 25 27 35 36 38 46 47 48 58 67 78 49/8 12 13 17 25 26 35 37 45 46 48 58 67 68 78 -77/4 12 14 18 23 28 36 37 46 47 56 57 58 68 78 49/8 12 14 17 23 27 35 38 46 48 57 58 67 68 78 -49/8 12 13 15 26 27 35 36 45 47 48 58 67 68 78 -7 12 13 15 26 28 35 37 45 46 47 58 67 68 78 -7/4 12 13 18 24 28 35 38 46 47 57 58 67 68 78 7 12 14 18 23 26 36 38 47 48 56 57 58 67 78 -7"





fullCycle = [[0,[]]]
i=0
graphInd = 0
edgeCount = 0
while i < len(inp):
    if edgeCount < 14:
        e1 = int(inp[i])-1
        e2 = int(inp[i+1])-1
        fullCycle[graphInd][1] += [[e1,e2]]
        i+=3
        edgeCount +=1
    else:
        j = i+1
        while j < len(inp) and inp[j] != " ":
            j+=1
        coef = sp.simplify(inp[i:j])
        fullCycle[graphInd][0] = coef
        i = j+1
        graphInd+=1
        fullCycle += [[0,[]]]
        edgeCount = 0

fullCycle = fullCycle[:-1]
print(fullCycle)

C = []
for G in fullCycle:
    nxG = readGraphfromEdgeList(G[1])
    C += [[G[0],nxG]]

pos8Cycle = lambda x: nx.circular_layout(nx.cycle_graph(8),center=x, scale=100)
pltChain(C,pos8Cycle,path=path+"sevenWheel")

for coef,G in C:
    if sum(np.array(G.degree)[:,1] != 3) > 2:
        C.remove([coef,G])

pltChain(C,pos8Cycle,path=path+"sevenWheelReduced")

gamma2 = psiWillwacher(C)
gamma2 = resolveMarkedIsos(gamma2)
gamma2 = resolveMarkedVanishing(gamma2)
pltChain(gamma2,pos8Cycle,path=path+"sevenWheelGamma2")

T = forgetNonTreePart(gamma2)
pltChain(T,pos8Cycle,path=path+"sevenWheelT")

tder = orientGraphForSder(T)
tder = resolveDirectedMarkedIsos(tder)
tder = resolveDirectedMarkedVanishing(tder)
pltChain(tder,pos8Cycle,path=path+"sevenWheelTder")

phiElement = step6Willwacher(T)
print(phiElement)