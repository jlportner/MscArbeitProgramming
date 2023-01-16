import numpy as np
import networkx as nx
from networkx.algorithms import isomorphism as nxiso
from sympy.combinatorics import Permutation

def getOrder(G):
    o = nx.get_edge_attributes(G, 'order')
    return {key: val for key, val in o.items() if val != -1}

def getEdgePerm(F1, F2, g):
    F1 = {(k[0], k[1]): v for k, v in F1.items()}
    F2 = {(k[0], k[1]): v for k, v in F2.items()}
    sig = np.zeros(len(F1))
    for e, i in F1.items():
        se = (g[e[0]], g[e[1]])
        if se in F2:
            sig[i - 1] = F2[se]
        else:
            sig[i - 1] = F2[(se[1], se[0])]
    perm = Permutation(sig - 1.0)
    return (-1) ** perm.parity()

"""
def getNodePerm(g,attributes):
    g = np.array(list(g.items()))
    sortInd = g[:,0].argsort()
    gSorted = g[sortInd]
    perm = Permutation(gSorted[:,1])
    return (-1) ** perm.parity()
"""

def resolveIsos(dG, withCof = True, node_match=None,edge_match=None, graphMatcher=nxiso.MultiGraphMatcher):
    i = 0
    dG = dG.copy()
    n = len(dG)
    while i < n - 1:
        H1 = dG[i]
        j = i + 1
        while j < n:
            H2 = dG[j]
            GM = graphMatcher(H1[1], H2[1],node_match = node_match,edge_match=edge_match)
            if GM.is_isomorphic():
                sgn = getEdgePerm(getOrder(H1[1]), getOrder(H2[1]), GM.mapping)
                #sgn *= getNodePerm(GM.mapping)
                H1[0] += sgn * H2[0]
                del dG[j]
                n += -1
            else:
                j += 1
        i += 1

    # remove 0s

    if withCof:
        result = [G for G in dG if G[0] != 0]
    else:
        result = dG
    return result

def resolveVanishing(dG,node_match=None,edge_match=None,graphMatcher=nxiso.MultiGraphMatcher):
    def hasOddAuto(H):
        GM = graphMatcher(H, H,node_match = node_match,edge_match=edge_match)
        for g in GM.isomorphisms_iter():
            F = getOrder(H)
            if getEdgePerm(F, F, g) == -1:
                return 1
        return 0

    result = [G for G in dG if not hasOddAuto(G[1])]
    return result