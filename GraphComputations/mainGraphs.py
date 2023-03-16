import networkx as nx
import numpy as np
from GCboundary import *


def createGamma2dash():
    H = nx.cycle_graph(6,nx.MultiGraph)
    H.add_edges_from([[1,5],[1,4],[2,4],[0,3]])
    for i,e in enumerate(H.edges):
        H.edges[e]['order'] = i+1

    G = nx.wheel_graph(6,nx.MultiGraph)
    for i,e in enumerate(G.edges):
        G.edges[e]['order'] = i+1

    C = [[2,G],[-5,H]]

    C = mark1Vertex(C)
    C = resolveMarkedIsos(C)
    C = resolveMarkedVanishing(C)
    C = split1MarkedIn2(C)
    C = resolveMarkedIsos(C)
    C = resolveMarkedVanishing(C)

    #dC = deltaMarkedGraphs(C)
    posFunc = lambda x: nx.circular_layout(nx.cycle_graph(9),center=x, scale=100)
    pltChain(C,posFunc,path="gamma2dash")
    #pltChain(dC, posFunc)
    np.save("g2d5wheel",C,allow_pickle=True)

#"""
createGamma2dash()

Cold = np.load("g2d5wheel.npy",allow_pickle=True).tolist()

H = nx.cycle_graph(6,nx.MultiGraph)
H.add_edges_from([[1,5],[1,4],[2,4],[0,3]])
for i,e in enumerate(H.edges):
    H.edges[e]['order'] = i+1

G = nx.wheel_graph(6,nx.MultiGraph)
for i,e in enumerate(G.edges):
    G.edges[e]['order'] = i+1


CO = [[2.0j,G],[-5.0j,H]]
C = psiWillwacher(CO)
C = resolveMarkedIsos(C)
C = resolveMarkedVanishing(C)
dC = deltaMarkedGraphs(C)
dC = resolveMarkedIsos(dC)
dC = resolveMarkedVanishing(dC)

Cold = np.load("g2d5wheel.npy",allow_pickle=True).tolist()

dM = Cold + dC
dM = resolveMarkedIsos(dM)
dM = resolveMarkedVanishing(dM)


posFunc = lambda x: nx.circular_layout(nx.cycle_graph(9),center=x, scale=100)
pltChain(dM,posFunc)
pltChain(C,posFunc,path="out2")
exit()

"""

G = nx.wheel_graph(4,nx.MultiGraph)
for i,e in enumerate(G.edges):
    G.edges[e]['order'] = i+1

CO = [[1,G]]
C = psiWillwacher(CO)
C = resolveMarkedIsos(C)
C = resolveMarkedVanishing(C)
dC = deltaMarkedGraphs(C)
dC = resolveMarkedIsos(dC)
dC = resolveMarkedVanishing(dC)

C = mark1Vertex([[1.0j,G]])
C = resolveMarkedIsos(C)
C = resolveMarkedVanishing(C)
C = split1MarkedIn2(C)
C = resolveMarkedIsos(C)
C = resolveMarkedVanishing(C)

dM = dC + C
dM = resolveMarkedIsos(dM)
posFunc = lambda x: nx.circular_layout(nx.cycle_graph(6),center=x, scale=100)
pltChain(dM,posFunc)
"""

"""
G = nx.wheel_graph(4,nx.MultiGraph)
for i,e in enumerate(G.edges):
    G.edges[e]['order'] = i+1

G1 = mark1Vertex([[1,G]])
G2d = split1MarkedIn2(G1)
G2 = psiWillwacher([[1,G]])
G2 = resolveMarkedIsos(G2)
G2 = resolveMarkedVanishing(G2)
posFunc = lambda x: nx.circular_layout(nx.cycle_graph(6),center=x, scale=100)
pltChain([[1,G]] + G1 + G2d + G2,posFunc)
"""