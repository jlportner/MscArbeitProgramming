import itertools

import networkx as nx
import numpy as np
from operator import itemgetter
from isoResolving import *
import sympy as sp

#Or do that they have to be -1 TODO
def nodeMatchMarked(u,v):
    return (( "extNr" not in u) and ("extNr" not in v)) or (( "extNr" in u) and ("extNr" in v) and u["extNr"] == v["extNr"])

def resolveMarkedIsos(dG, withCof = True):
    return resolveIsos(dG, withCof, node_match=nodeMatchMarked)

def resolveMarkedVanishing(dG):
    return resolveVanishing(dG,node_match=nodeMatchMarked)

def getMarkedNodes(G):
    o = nx.get_node_attributes(G, 'extNr')
    return {key: val for key, val in o.items() if val != -1}

def mark1Vertex(C):
    dC = []
    for cof, G in C:
        dG = []
        #permFactor = np.math.factorial(len(G.nodes)-1)
        for v in G.nodes:
            H = G.copy()
            #H.nodes[v]["ext"] = True
            H.nodes[v]["extNr"] = 1
            dG += [[cof,H]]
        dG = resolveMarkedIsos(dG)
        GM = nxiso.GraphMatcher(G, G, node_match=nodeMatchMarked)
        #autSize = len(list(GM.isomorphisms_iter()))
        #dG = [[cof/autSize,G] for cof,G in dG]
        dC += dG

    return dC

#assume ordering on G and H starts at 1, only works on simple graphs
def graphComposition(G,H,j,coef=1):
    if j not in G.nodes:
        raise ValueError("Node not in graph")

    G = G.copy()
    H = H.copy()

    G.nodes[j]["speziale"] = True

    #adapt ordering on H to be after G
    orderOffset = len(G.edges)
    for e in H.edges:
        H.edges[e]['order'] += orderOffset

    Q = nx.disjoint_union(G,H)
    jNew = list(nx.get_node_attributes(Q,"speziale").keys())[0]

    nbhd = np.array(list(Q.edges(jNew, keys=True)))
    nbhSize = len(nbhd)
    Hlabels = np.arange(len(G),len(G)+len(H))
    dG = []
    for perm in itertools.product(Hlabels, repeat=nbhSize):
        newEdges = nbhd.copy()
        newEdges[:, 0] = perm
        newQ = Q.copy()
        newQ.remove_node(jNew)
        for i,e in enumerate(newEdges):
            newQ.add_edge(*e)
            data = Q.get_edge_data(*nbhd[i])
            newQ.edges[e].update(data)

        newQ = nx.convert_node_labels_to_integers(newQ)
        dG += [[coef,newQ]]

    return dG

def split1MarkedIn2(C):

    H = nx.MultiGraph()
    H.add_nodes_from([0,1])
    H.nodes[0]["extNr"] = 1
    H.nodes[1]["extNr"] = 2

    dG = []
    for cof, G in C:
        markedNodes = list(nx.get_node_attributes(G, "extNr").keys())
        if len(markedNodes) != 1:
            raise ValueError("Graph needs to have excatly one marked vertex")

        dG += graphComposition(G,H,markedNodes[0],coef=cof)

    dG = resolveMarkedIsos(dG)
    dG = resolveMarkedVanishing(dG)

    dC = []
    for cof, Q in dG:
        markedV = set(nx.get_node_attributes(Q, "extNr").keys())
        unmarkedV = set(Q.nodes) - markedV
        degrees = dict(Q.degree)
        if np.min(itemgetter(*markedV)(degrees)) <= 0 or np.min(itemgetter(*unmarkedV)(degrees)) < 3:
            continue
        else:
            dC += [[cof, Q]]

    return dC

def deltaMarkedGraphs(C):
    HExt = nx.path_graph(2,nx.MultiGraph)
    HExt.nodes[0]["extNr"] = 1
    HExt.edges[(0,1,0)]["order"] = 1

    HInt = nx.path_graph(2,nx.MultiGraph)
    HInt.edges[(0, 1,0)]["order"] = 1

    dG = []
    for cof, G in C:
        for v in G.nodes:
            if "extNr" in G.nodes[v]:
                if G.degree[v] >= 2:
                    HExt.nodes[0]["extNr"] = G.nodes[v]["extNr"]
                    out= graphComposition(G,HExt,v,coef=cof)
                    #for el in out:
                    #    el[0] *= 1
                    dG += out
            else:
                if G.degree[v] >= 4:
                    out = graphComposition(G, HInt, v, coef=cof)
                    for el in out:
                        el[0] *= 1/2
                    dG += out

    dC = []
    for cof, H in dG:
        markedV = set(nx.get_node_attributes(H, "extNr").keys())
        unmarkedV = set(H.nodes) - markedV
        degrees = dict(H.degree)
        if np.min(itemgetter(*unmarkedV)(degrees)) < 3:
            continue
        else:
            dC += [[cof, H]]

    return dC

def psiWillwacher(C):
    dC = []
    for cof, G in C:
        #GM = nxiso.GraphMatcher(G, G, node_match=nodeMatchMarked)
        #autSize = len(list(GM.isomorphisms_iter()))

        for e in G.edges:
            H = G.copy() #type: nx.MultiGraph
            orderE = H.edges[e]['order']
            ord = getOrder(H)
            for e2 in H.edges:
                if H.edges[e2]['order'] > orderE:
                    H.edges[e2]['order'] -=1
            H.remove_edge(*e)
            H2 = H.copy()
            H.nodes[e[0]]['extNr'] = 1
            H.nodes[e[1]]['extNr'] = 2

            H2.nodes[e[0]]['extNr'] = 2
            H2.nodes[e[1]]['extNr'] = 1

            #dC += [[sp.Rational(cof)/autSize * (-1)**(orderE-1),H],[sp.Rational(cof)/autSize * (-1)**(orderE-1),H2]]
            dC += [[cof * (-1)**(orderE-1),H],[cof * (-1)**(orderE-1),H2]]

    return dC

def forgetNonTreePart(C):
    def isInternal3Tree(G: nx.MultiGraph) -> bool:
        H = G.copy()
        extNodes = set(nx.get_node_attributes(H, "extNr").keys())
        intNodes = set(H.nodes) - extNodes
        internalDegrees = itemgetter(*intNodes)(dict(H.degree))

        if np.any(internalDegrees != np.array(3)):
            return False
        H.remove_nodes_from(list(extNodes))
        result = nx.is_tree(H)
        return result

    result = [G for G in C if isInternal3Tree(G[1])]
    return result
