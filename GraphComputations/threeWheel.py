import networkx as nx
import numpy as np
from GCboundary import *
from tderIdentification import *

G = nx.wheel_graph(12,nx.MultiGraph)
for i,e in enumerate(G.edges):
    G.edges[e]['order'] = i+1


C = [[1,G]]

#pos6Cycle = lambda x: nx.circular_layout(nx.cycle_graph(6),center=x, scale=100)

gamma2 = psiWillwacher(C)
gamma2 = resolveMarkedIsos(gamma2)
gamma2 = resolveMarkedVanishing(gamma2)

T = forgetNonTreePart(gamma2)

tder = orientGraphForSder(T)
tder = resolveDirectedMarkedIsos(tder)
tder = resolveDirectedMarkedVanishing(tder)

phiElement = step6Willwacher(T)
print(phiElement)