import networkx as nx
import numpy as np
from GCboundary import *
from tderIdentification import *

H = nx.cycle_graph(6,nx.MultiGraph)
H.add_edges_from([[1,5],[1,4],[2,4],[0,3]])
for i,e in enumerate(H.edges):
    H.edges[e]['order'] = i+1

G = nx.wheel_graph(6,nx.MultiGraph)
for i,e in enumerate(G.edges):
    G.edges[e]['order'] = i+1

C = [[-5,H]]

gamma2 = psiWillwacher(C)
gamma2 = resolveMarkedIsos(gamma2)
gamma2 = resolveMarkedVanishing(gamma2)
T = forgetNonTreePart(gamma2)
phiElement = step6Willwacher(T)
print(phiElement)