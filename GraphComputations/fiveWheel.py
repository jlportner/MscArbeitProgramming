import networkx as nx
import numpy as np
from GCboundary import *
from sage.all import *
from tderIdentification import *
path = "/media/Transfer/ETH/MscArbeit/Images/"

H = nx.cycle_graph(6,nx.MultiGraph)
H.add_edges_from([[1,5],[1,4],[2,4],[0,3]])
for i,e in enumerate(H.edges):
    H.edges[e]['order'] = i+1

G = nx.wheel_graph(6,nx.MultiGraph)
for i,e in enumerate(G.edges):
    G.edges[e]['order'] = i+1

nx.write_edgelist(G,"fiveWheel")

C = [[1,G],[-Rational(5/2),H]]

pos6Cycle = lambda x: nx.circular_layout(nx.cycle_graph(6),center=x, scale=100)
pos7Cycle = lambda x: nx.circular_layout(nx.cycle_graph(7),center=x, scale=100)
pltChain(C,pos6Cycle,path=path+"fiveWheel")

gamma1 = mark1Vertex(C)
gamma1 = resolveMarkedIsos(gamma1)
gamma1 = resolveMarkedVanishing(gamma1)
pltChain(gamma1,pos6Cycle,path=path+"fiveWheelGamma1")

gamma2Dash = split1MarkedIn2(gamma1)
gamma2Dash = resolveMarkedIsos(gamma2Dash)
gamma2Dash = resolveMarkedVanishing(gamma2Dash)
pltChain(gamma2Dash,pos7Cycle,path=path+"fiveWheelGamma2Dash")

gamma2 = psiWillwacher(C)
gamma2 = resolveMarkedIsos(gamma2)
gamma2 = resolveMarkedVanishing(gamma2)
pltChain(gamma2,pos6Cycle,path=path+"fiveWheelGamma2")

T = forgetNonTreePart(gamma2)
pltChain(T,pos6Cycle,path=path+"fiveWheelT")

tder = orientGraphForSder(T)
tder = resolveDirectedMarkedIsos(tder)
tder = resolveDirectedMarkedVanishing(tder)
pltChain(tder,pos6Cycle,path=path+"fiveWheelTder")

phiElement = step6Willwacher(T)
print(phiElement)

res = deltaMarkedGraphs(gamma2)
res = resolveMarkedIsos(res)
res = resolveMarkedVanishing(res)

#res = [[-cof,G] for (cof,G) in res]
res += gamma2Dash
res = resolveMarkedIsos(res)
print(res)
#pltChain(res,pos7Cycle,path="out")