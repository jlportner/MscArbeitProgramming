import networkx as nx
import numpy as np
from GCboundary import *
from tderIdentification import *
path = "/media/Transfer/ETH/MscArbeit/Images/"

def wheelLayout(n,center = (0,0)):
    pos = {0 : np.array(center)}
    scale = 100.0
    lAng = np.linspace(0, 2*np.pi, n+1)[:-1]
    lPoints = np.vstack((scale * np.cos(lAng) + center[0], scale * np.sin(lAng) + center[1])).T
    pos |= dict(enumerate(lPoints.tolist(), 1))
    return pos


G = nx.wheel_graph(4,nx.MultiGraph)
for i,e in enumerate(G.edges):
    G.edges[e]['order'] = i+1

#H = nx.wheel_graph(14,nx.MultiGraph)
#for i,e in enumerate(G.edges):
#    G.edges[e]['order'] = i+1

C = [[1,G]]

posWheel = lambda x: wheelLayout(3,center=(600,0))
pltChain(C,posWheel,path="out")

#exit()

gamma2 = psiWillwacher(C)
gamma2 = resolveMarkedIsos(gamma2)
gamma2 = resolveMarkedVanishing(gamma2)
pltChain(gamma2,posWheel,path="out2")
exit()

T = forgetNonTreePart(gamma2)

tder = orientGraphForSder(T)
tder = resolveDirectedMarkedIsos(tder)
tder = resolveDirectedMarkedVanishing(tder)

phiElement = step6Willwacher(T)
print(phiElement)