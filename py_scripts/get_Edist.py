import networkx as nx
import numpy as np
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

nnodes = H.number_of_nodes()

#define two empty distance matrices
# 1 is the energy difference between source min and highest TS on path
# 2 is the sum of all barriers 
distances = np.zeros((nnodes,nnodes))
distances2 = np.zeros((nnodes,nnodes))

other_dists = dict()
other_dists2 = dict()

for node in H.nodes():
    other_dists[node] = dict()
    other_dists2[node] = dict()    

# the node id is not identical to the array id, so we kepe track of it 
i = -1
#total number of distances calculated
total_dists = 0

for inode in H.nodes():
    i += 1
    j = -1
    for jnode in H.nodes():
        # initialise everything here
        j += 1
        E = list()
        E2 = list()
        calculated_dists = 0
        # if we go from i to i, the distance is 0
        if (i==j):
            distances[i][j] = 0.0
            distances2[i][j] = 0.0
            calculated_dists += 1
        # if the key jnode exists in the distance dictionary for inode, then we already calculated this distance, no need to redo that calculation 
        elif jnode in other_dists[inode]:
            distances[i][j] = other_dists[inode][jnode]
            distances2[i][j] = other_dists2[inode][jnode]
        # otherwise we need to find the shortest path
        else:
            # get the dijkstra path using barrier height as the weight
            path = nx.dijkstra_path(H,inode,jnode)
            nedges = len(path) - 1
            # iterate over all pairs of adjacent ndoes to yield the edges, and store the energies we are interested in
            for k in range(nedges):
                min1 = path[k]
                min2 = path[k+1]
                E.append(H[min1][min2]["E"])        # the absolute E of the TS
                E2.append(H[min1][min2]["weight"])  # the barrier height
            # at this point we have list of energies along the path
            for k in range(nedges):
                for l in range(k+1,nedges+1):
                    #k is the first index, l the last, E[k:l] gives the barriers between them
                    min1 = path[k]
                    min2 = path[l]
                    # check whether we already hav that distance (we might not have the overall path, but subsections of it)
                    if not(min2 in other_dists[min1]):
                        refE = H.nodes()[min1]['E']
                        E1 = max(E[k:l]) - refE
                        E2sum = np.sum(np.asarray(E2[k:l]))
                        other_dists[min1][min2] = E1
                        other_dists2[min1][min2] = E2sum
                        if (min1==inode) and (min2==jnode):
                            distances[i][j] = E1
                            distances2[i][j] = E2sum
                        calculated_dists += 1
        total_dists += calculated_dists
        print("Completed shortest path calculation for minima ", i+1, " and ", j+1)
        print("Calculated ", calculated_dists, " new distances this iteration, total distances found: ", total_dists)

np.savetxt("dist1.csv",distances,delimiter = ',')
np.savetxt("dist2.csv",distances2,delimiter = ',')

outf = open("min.retained", "w")
outf.write(str(len(H.nodes())) +  "\n")
for node in H.nodes():
    outf.write("{0:7d} \n".format(node))