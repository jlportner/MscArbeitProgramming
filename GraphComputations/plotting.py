import os
import networkx as nx
from isoResolving import *
import matplotlib.pyplot as plt
from markedGraphs import *

def pltGraph(H, pos,path = "out"):
    G = H.copy()
    for n, npos in pos.items():
        if(G.has_node(n)):
            G.nodes[n]['pos'] = '"%d,%d"'%(npos[0],npos[1])
            G.nodes[n]['shape'] = "point"
    p = nx.drawing.nx_pydot.to_pydot(G)
    p.write(path + ".dot")
    os.system("neato -n2 -Tpng " + path + ".dot -o " + path + ".png")

def pltChain(C,posFunction,path="out",lineWidth=5):
    if nx.is_directed(C[0][1]):
        masterG = nx.MultiDiGraph()
    else:
        masterG = nx.MultiGraph()
    hPos = 0
    for j, (count,H) in enumerate(C):
        if j!= 0 and j % lineWidth == 0:
            hPos -= 3.0
        vPos = 3 * (j % lineWidth)
        pos = posFunction((100*vPos, 100*hPos))
        G = H.copy()
        nx.set_edge_attributes(G,2.0,name="penwidth")
        edgeOrder = getOrder(G)
        nx.set_edge_attributes(G, {key: {"label": val} for key, val in edgeOrder.items()})

        for node, nodePos in pos.items():
            if(G.has_node(node)):
                G.nodes[node]['pos'] = '"%d,%d"' % (nodePos[0], nodePos[1])
                G.nodes[node]['shape'] = "point"
                G.nodes[node]['width'] = "0.15pt"

        markedNodes = getMarkedNodes(G)
        nx.set_node_attributes(G, {key: {"label": val, "fontsize": "12pt", "shape": "circle", "fixedsize": True, "width": "0.25pt", "penwidth": "2pt", "color": "#ff0000"} for
                                   key, val in markedNodes.items()})

        masterG = nx.disjoint_union(masterG,G)
        cof = str(j) + "coef"
        masterG.add_node(cof)
        if j!= 0 and not (j % lineWidth == 0):
            if count >= 0:
                masterG.nodes[cof]['label'] = "+ " + str(count)
            else:
                masterG.nodes[cof]['label'] = "− " + str(abs(count))
        else:
            masterG.nodes[cof]['label'] = count
        masterG.nodes[cof]['shape'] = "plaintext"
        masterG.nodes[cof]['fontsize'] = "24pt"
        masterG.nodes[cof]['pos'] = '"%d,%d"' %(100*(-1.45 + vPos), hPos*100)

    nx.drawing.nx_pydot.write_dot(masterG,path + ".dot")

    os.system("neato -n2 -Tpdf " + path + ".dot -o " + path + ".pdf")

def pltStandardGraph(C,posFunction,path="out",lineWidth=5):
    masterG = nx.MultiGraph()
    hPos = 0
    for j, H in enumerate(C):
        if j!= 0 and j % lineWidth == 0:
            hPos -= 3.0
        vPos = 6 * (j % lineWidth)
        pos = posFunction((vPos, hPos))
        G = H.copy()
        nx.set_edge_attributes(G,2.0,name="penwidth")

        for node, nodePos in pos.items():
            if(G.has_node(node)):
                G.nodes[node]['pos'] = '"%d,%d"' % (nodePos[0], nodePos[1])
                G.nodes[node]['shape'] = "point"
                G.nodes[node]['width'] = "0.15pt"
        masterG = nx.disjoint_union(masterG,G)


    p = nx.drawing.nx_pydot.to_pydot(masterG)

    p.write(path + ".dot")
    os.system("neato -n2 -Tpdf " + path + ".dot -o " + path + ".pdf")