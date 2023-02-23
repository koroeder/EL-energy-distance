import pickle
import networkx as nx
import numpy as np
import os
import json

# change str keys to int keys
def str2int_dics(d):
    out = dict()
    for k1,v1 in d.items():
        out[int(k1)] = dict()
        for k2, v2 in v1.items():
            out[int(k1)][int(k2)] = v2
    return out

#open pickled graph
with open('ELgraph.gpickle', 'rb') as f:
    G = pickle.load(f)

H = nx.subgraph(G, list(sorted(nx.strongly_connected_components(G), key=len, reverse=True))[0])
print(H)
nodes = H.nodes()
nnodes = H.number_of_nodes()

#define two empty distance matrices
# 1 is the energy difference between source min and highest TS on path
# 2 is the sum of all barriers 
distances = np.zeros((nnodes,nnodes))
distances2 = np.zeros((nnodes,nnodes))

other_dists = dict()
other_dists2 = dict()

for node in nodes:
    other_dists[node] = dict()
    other_dists2[node] = dict()    

# check whether we have a backup directory, if we don't create it
bupath = "./backup/"
if not os.path.exists(bupath):
    os.makedirs(bupath)
    idx = -1
# if we do, we check what the index of the last calculated node is
else:
    for idx in range(nnodes):
        if os.path.exists(bupath+"dist1."+str(idx+1)+".csv"):
            distances[idx][:] = np.loadtxt(bupath+"dist1."+str(idx+1)+".csv", delimiter = ',')
            distances2[idx][:] = np.loadtxt(bupath+"dist2."+str(idx+1)+".csv", delimiter = ',')
        else:
            print("Opening backup files up to index ", idx+1)
            with open("bu_dist1.pkl", "r") as f:
                other_dists = json.load(f)
            with open("bu_dist2.pkl", "r") as f:
                other_dists2 = json.load(f) 
            other_dists = str2int_dics(other_dists)
            other_dists2 = str2int_dics(other_dists2)
            break

#total number of distances calculated
total_dists = 0
# the node id is not identical to the array id, so we kepe track of it 
i = -1
for inode in nodes:    
    i += 1
    if i>=idx:
        j = -1
        for jnode in nodes:
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
            if ((j+1)%1000==0):
                print("Completed shortest path calculation for minima ", i+1, " and ", j+1)
                print("Calculated ", calculated_dists, " new distances this iteration, total distances found: ", total_dists)
        np.savetxt(bupath+"dist1."+str(i+1)+".csv", distances[i][:], delimiter = ',')
        np.savetxt(bupath+"dist2."+str(i+1)+".csv", distances2[i][:], delimiter = ',')
        
        with open("bu_dist1.pkl", "w") as f:
            json.dump(other_dists, f)
        with open("bu_dist2.pkl", "w") as f:
            json.dump(other_dists2, f)     
    else:
        print("Skipping calculations for index ", i)   

np.savetxt("dist1.csv", distances, delimiter = ',')
np.savetxt("dist2.csv", distances2, delimiter = ',')

