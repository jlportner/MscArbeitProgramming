import networkx as nx
import numpy as np
import more_itertools as mit
from plotting import *
from isoResolving import *

def deltaGC(C):
    dC = []
    for cof,G in C:
        dG = []
        newNodeLabel = min(set(np.arange(max(G.nodes)+2)) - set(G.nodes))
        for v in G.nodes:
            nbhd = list(G.edges(v,keys=True))
            nbhSize = len(nbhd)

            if nbhSize > 3:
                for part in mit.set_partitions(nbhd,2):
                    if len(part[0]) == nbhSize or len(part[0]) < 2 or len(part[1]) < 2:
                        continue
                    for j in range(2):
                        P = [part[j],part[j-1]]
                        H = G.copy() # type: nx.MultiGraph
                        H.remove_node(v)
                        u,w = v, newNodeLabel #str(v) + 'a', str(v) + 'b' #TODO: The second argument needs fixing
                        H.add_nodes_from((u,w))
                        newNodes = [u,w]

                        for i in range(2):
                            for e in P[i]:
                                data = G.get_edge_data(*e)
                                if(e[0] != v):
                                    key = H.add_edge(newNodes[i],e[0])
                                    eNew = (newNodes[i],e[0],key)
                                else:
                                    key = H.add_edge(newNodes[i], e[1])
                                    eNew = (newNodes[i],e[1],key)
                                H.edges[eNew].update(data)

                        e = (u,w,0)
                        H.add_edge(*e)
                        H.edges[e]['order'] = len(G.edges) + 1
                        dG.append([cof,H])
        dC += dG
    return dC

