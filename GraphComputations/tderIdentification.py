import subprocess
import networkx as nx
import numpy as np
from GCboundary import *


def resolveDirectedMarkedIsos(dG, withCof = True):
    return resolveIsos(dG, withCof, node_match=nodeMatchMarked,graphMatcher=nxiso.MultiDiGraphMatcher)

def resolveDirectedMarkedVanishing(dG):
    return resolveVanishing(dG,node_match=nodeMatchMarked,graphMatcher=nxiso.MultiDiGraphMatcher)

def orientRootedTree(G,root,inComing,minLabel):
    """

    :type G: nx.MultiDiGraph
    """
    newLabel = minLabel

    if "extNr" in G.nodes[root]:
        return 0

    incomingEdges = [ed for ed in G.in_edges if (ed[1] == root)]
    for e in incomingEdges:
        if e[0] == inComing:
            continue
        G.remove_edge(*e)
        G.edges[(e[1],e[0],e[2])]["order"] = newLabel
        labelIncrease = orientRootedTree(G,e[0],root,newLabel+1)
        newLabel += labelIncrease + 1


    return newLabel - minLabel

def orientGraphForSder(C):
    result = []
    for cof, G in C: # type: int, nx.MultiGraph
        extVert = getMarkedNodes(G)
        for v in extVert:
            nbhdEdges = [ed for ed in G.edges if (ed[0] == v or ed[1] == v)]
            for e in nbhdEdges:
                dirG = G.to_directed() #type: nx.MultiDiGraph
                u = e[0]
                if e[0] == v:
                    u = e[1]
                dirG.remove_edge(u,v,e[2])
                dirG.edges[(v,u,e[2])]["order"] = 1
                orientRootedTree(dirG,u,v,2)

                H = dirG.to_undirected()

                sig = np.zeros(len(G.edges))
                for e,i in getOrder(G).items():
                    sig[i-1] = H.edges[e]["order"]
                perm = Permutation(sig - 1.0)
                sgn = (-1) ** perm.parity()
                result += [[cof * sgn,dirG]]

    return result

def treeToLieRecursion(G,root) -> str:
    """

        :type G: nx.MultiDiGraph
    """
    if "extNr" in G.nodes[root]:
        return "X" + str(G.nodes[root]["extNr"])

    outEdges = [ed for ed in G.in_edges if (ed[0] == root)]
    assert len(outEdges) == 2
    if G.edges[outEdges[0]]["order"] < G.edges[outEdges[1]]["order"]:
        w1,w2 = outEdges[0][1],outEdges[1][1]
    else:
        w1, w2 = outEdges[1][1], outEdges[0][1]

    w1Part = treeToLieRecursion(G, w1)
    w2Part = treeToLieRecursion(G, w2)
    return "L([" + w1Part + "," + w2Part + "])"

def tderToLie(C):
    k = len(getMarkedNodes(C[0][1]))
    result = [None] * k

    for cof, G in C:  # type: int, nx.MultiDiGraph
        extVert = getMarkedNodes(G)
        extOutEdges = [ed for ed in G.edges if (ed[0] in extVert.keys())]
        assert len(extOutEdges) == 1
        root = extOutEdges[0][1]
        origin = extOutEdges[0][0]
        index = G.nodes[origin]["extNr"]
        expr = treeToLieRecursion(G,root)
        if result[index-1] is None:
            result[index - 1] = str(cof) + "*" + expr
        else:
            result[index-1] += "+" + str(cof) + "*" + expr

    return result

def step6Willwacher(C):
    C = orientGraphForSder(C)
    C = resolveDirectedMarkedIsos(C)
    C = resolveDirectedMarkedVanishing(C)
    tderEl = tderToLie(C)
    assert len(tderEl) == 2
    expr = tderEl[0].replace("X2","Y").replace("X1","X") + " - (" + tderEl[1].replace("X2","Y").replace("X1","X") + ")"
    return simplifyLieExpressionWithSage(expr,["Y","X"])

def simplifyLieExpressionWithSage(expr, vars):
    cmd = "L.<"
    for x in vars[:-1]:
        cmd += str(x) + ","

    cmd += str(vars[-1]) + "> = LieAlgebra(QQ,abelian=False); print(" + expr + ")"

    result = subprocess.run(["sage","-c",cmd], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    simpleExpr = result.stdout
    return str(simpleExpr)[2:-3]
