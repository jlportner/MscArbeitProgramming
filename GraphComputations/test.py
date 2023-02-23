import networkx as nx
import numpy as np
from GCboundary import *
from tderIdentification import *

path = "/media/Transfer/ETH/MscProgramming/py/outFiles/"
H = nx.cycle_graph(6,nx.MultiGraph)
H.add_edges_from([[1,5],[1,4],[2,4],[0,3]])
for i,e in enumerate(H.edges):
    H.edges[e]['order'] = i+1

G = nx.wheel_graph(4,nx.MultiGraph)
for i,e in enumerate(G.edges):
    G.edges[e]['order'] = i+1

#C = [[-5,H]]
C = [[1,G]]
#C = [[2,G],[-5,H]]
"""
pos6Cycle = lambda x: nx.circular_layout(nx.cycle_graph(6),center=x, scale=100)
gamma2 = psiWillwacher(C)
gamma2 = resolveMarkedIsos(gamma2)
gamma2 = resolveMarkedVanishing(gamma2)
pltChain(gamma2,pos6Cycle,path=path+"GGamma2")

T = forgetNonTreePart(gamma2)
pltChain(T,pos6Cycle,path=path+"GT")

tder = orientGraphForSder(T)
tder = resolveDirectedMarkedIsos(tder)
tder = resolveDirectedMarkedVanishing(tder)
pltChain(tder,pos6Cycle,path=path+"GTder")
"""
gamma2 = psiWillwacher(C)
gamma2 = resolveMarkedIsos(gamma2)
gamma2 = resolveMarkedVanishing(gamma2)
T = forgetNonTreePart(gamma2)
phiElement = step6Willwacher(T)
print(phiElement)