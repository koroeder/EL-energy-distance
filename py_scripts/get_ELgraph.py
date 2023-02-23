import networkx as nx
import numpy as np
import pickle
import sys
import math

#create graph from min.data and ts.data
mins = np.genfromtxt("min.data", dtype = float, usecols=(0,))
nmin = len(mins)
G = nx.DiGraph()

for vert in range(nmin):
    G.add_node(vert+1, E=mins[vert])

tiny = 0.000001

with open("ts.data","r") as f:
    for line in f:
        Ets = float(line.split()[0])
        min1 = int(line.split()[3])
        E1 = mins[min1-1]
        min2 = int(line.split()[4])
        E2 = mins[min2-1]
        if (Ets<E1) or (Ets<E2):
            Ets = max(E1,E2) + tiny
        if (min1!=min2):
            if G.has_edge(min1,min2):
               if (G.edges[min1,min2]['E']>Ets):
                   G.edges[min1,min2]['E'] = Ets
                   G.edges[min1,min2]['weight'] = Ets-E1
                   G.edges[min2,min1]['E'] = Ets
                   G.edges[min2,min1]['weight'] = Ets-E2
            else:
                G.add_edge(min1, min2, weight=Ets-E1, E=Ets)
                G.add_edge(min2, min1, weight=Ets-E2, E=Ets)

print( 'Number of minima: ',G.number_of_nodes())
print( 'Number of transition states: ', int(G.number_of_edges()/2))

H = nx.subgraph(G, list(sorted(nx.strongly_connected_components(G), key=len, reverse=True))[0])

print(H)

outf = open("min.retained", "w")
outf.write(str(len(H.nodes())) +  "\n")
for node in H.nodes():
    outf.write("{0:7d} \n".format(node))

with open('ELgraph.gpickle','wb') as f:
    pickle.dump(G, f, pickle.HIGHEST_PROTOCOL)
