import networkx as nx
import numpy as np
import itertools
from plotting import *
import sympy as sp

def recursionSimpleGraph(G,freeDeg,orderStart=1):
    l = len(freeDeg)
    if l <= 3:
        if l == 0:
            return [nx.MultiGraph(G)]
        if l == 1:
            del G
            return []
        elif l == 2:
            if sum(freeDeg.values()) != 2:
                del G
                return []
            else:
                G.add_edges_from([list(freeDeg.keys())])
                #G.edges[list(freeDeg.keys())]['order'] = orderStart
                return [nx.MultiGraph(G)]

        else:
            val = np.array(list(freeDeg.items()))
            val= val[val[:, 1].argsort()]

            if np.all(val[:, 1] == [2, 2, 2]):
                G.add_edges_from([(val[0, 0], val[2, 0]), (val[1, 0], val[2, 0]), (val[1, 0], val[0, 0])])
                #print("1")
                return [nx.MultiGraph(G)]

            elif np.all(val[:, 1] == [1, 1, 2]):
                G.add_edges_from([(val[0, 0], val[2, 0]), (val[1, 0], val[2, 0])])
                return [nx.MultiGraph(G)]

            else:
                del G
                return []

    v = min(freeDeg)
    C = []
    nodes = list(freeDeg.keys())
    edges = np.vstack([v*np.ones_like(nodes),nodes]).T[1:]

    for E in itertools.combinations(edges.tolist(), freeDeg[v]):
        newDeg = freeDeg.copy()
        newDeg[v] = 0
        H = G.copy()
        H.add_edges_from(E)
        #order = orderStart
        for e in E:
            newDeg[e[1]] -= 1
            #H.edges[e]['order'] = order
            #order += 1

        newDeg = {d:newDeg[d] for d in newDeg if newDeg[d] > 0}
        C += recursionSimpleGraph(H,newDeg)#,orderStart=order)

    print(len(C))
    return C

def recursionMultiGraph(G,freeDeg):
    l = len(freeDeg)
    if l <= 2:
        if l == 0:
            return [nx.MultiGraph(G)]
        if l == 1:
            del G
            return []
        else:
            val = list(freeDeg.values())
            if val[0] != val[1]:
                del G
                return []
            else:
                for i in range(val[0]):
                    G.add_edges_from([list(freeDeg.keys())])
                return [nx.MultiGraph(G)]

    v = min(freeDeg)
    C = []
    nodes = list(freeDeg.keys())
    edges = np.vstack([v*np.ones_like(nodes),nodes]).T[1:]

    for E in itertools.combinations_with_replacement(edges.tolist(), freeDeg[v]):
        newDeg = freeDeg.copy()
        newDeg[v] = 0
        H = G.copy()
        H.add_edges_from(E)
        for e in E:
            newDeg[e[1]] -= 1
            if newDeg[e[1]] < 0:
                continue
        if not min(newDeg.values()) < 0:
            newDeg = {d:newDeg[d] for d in newDeg if newDeg[d] > 0}
            C += recursionMultiGraph(H,newDeg)

    print(len(C))
    return C



def createAll3ValGraph(rank):
    v = 2 * rank - 2
    freeDeg = {d:3 for d in range(v)}
    return createAllGraphsWithDeg(freeDeg)

def createAllGraphsFromList(li):
    v = len(li)
    degVec = dict(np.vstack((np.arange(v),li)).T.tolist())
    return createAllGraphsWithDeg(degVec)
def createAllGraphsWithDeg(degVec):
    v = len(degVec)
    G = nx.empty_graph(v,create_using=nx.MultiGraph)
    G.nodes[0]['extNr'] = 1
    G.nodes[1]['extNr'] = 2
    C = recursionSimpleGraph(G,degVec)
    print(len(C))
    C = [G for G in C if nx.is_connected(G) and not (0,1,0) in G.edges]
    for G in C:
        for i, e in enumerate(G.edges):
            G.edges[e]['order'] = i + 1
    return C

C = []
degList = [[2,2,4,4,3,3],[2,2,5,3,3,3],[3,3,3,3,3,3],
            #[2,1,4,4,4,3],[2,1,5,4,3,3],[2,1,6,3,3,3],
            #[3,1,4,4,3,3],[3,1,5,3,3,3],[3,2,4,3,3,3],
            #[4,1,4,3,3,3],[4,2,3,3,3,3],[5,1,3,3,3,3]
           ]
for degV in degList:
    C += createAllGraphsFromList(degV)

C = [[1,G] for G in C]
C = resolveMarkedIsos(C,withCof=False)
C = resolveMarkedVanishing(C)

varlist = []
for i,c in enumerate(C):
    x = sp.symbols('x' + str(i))
    c[0] = x
    varlist.append(x)

dC = deltaMarkedGraphs(C)
Cold = np.load("g2dH.npy",allow_pickle=True).tolist()
dC += Cold
dC = resolveMarkedIsos(dC)
dC = resolveMarkedVanishing(dC)
B1 = np.array(dC, dtype=object)
eqs = B1[:,0].tolist()
sol = sp.linsolve(eqs, varlist)

print(sol)


posFunc = lambda x: nx.circular_layout(nx.cycle_graph(9),center=x, scale=100)
pltChain(C + dC,posFunc,path="allGraphs")
