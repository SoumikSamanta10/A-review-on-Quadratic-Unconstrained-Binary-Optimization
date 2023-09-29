# import structure for Q matrix
from collections import defaultdict
from math import sqrt

# Import networkx for graph tools
import networkx as nx

# Import dwave_networkx for d-wave graph tools/functions
import dwave_networkx as dnx

# Import matplotlib.pyplot to draw graphs on screen
import matplotlib
matplotlib.use("agg")    # must select backend before importing pyplot
import matplotlib.pyplot as plt

# Set the solver we're going to use
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

sampler = EmbeddingComposite(DWaveSampler())

# Create empty graph
G = nx.Graph()

# Add edges to graph - this also adds the nodes
G.add_weighted_edges_from([(1, 2, 1), (1, 3, 1), (1, 4, 2), (2, 3, 4), (2, 4, 3), (3, 4, 3)])


Q = defaultdict(int)
P = 8
chainstrength = 8
numruns = 10

num_of_towns = len(G.nodes)

#print(num_of_towns)

counter = 0


for t1 in range(num_of_towns):  #time loop
    #print("\n{} \n".format(t1))
    for v1, v2, d in G.edges(data=True):  #town loop
        #print("{} {} {}\t".format(v1, v2, d["weight"]))
        
        #diagonal terms
        Q[(num_of_towns * t1 + v1, num_of_towns * t1 + v1)] = -2*P
        Q[(num_of_towns * t1 + v2, num_of_towns * t1 + v2)] = -2*P

        # different vertex at same time
        Q[(num_of_towns * t1 + v1, num_of_towns * t1 + v2)] += P
        Q[(num_of_towns * t1 + v2, num_of_towns * t1 + v1)] += P

        # acceptable transitions
        Q[(num_of_towns * t1 + v1, num_of_towns * ((t1+1)%num_of_towns) + v2)] += d["weight"]
        Q[(num_of_towns * t1 + v2, num_of_towns * ((t1+1)%num_of_towns) + v1)] += d["weight"]

        #0 terms
        Q[(num_of_towns * t1 + v1, num_of_towns * ((t1+2)%4) + v2)] += 0 
        Q[(num_of_towns * t1 + v2, num_of_towns * ((t1+2)%4) + v1)] += 0
        Q[(num_of_towns * t1 + v1, num_of_towns * ((t1+3)%4) + v2)] += 0
        Q[(num_of_towns * t1 + v2, num_of_towns * ((t1+3)%4) + v1)] += 0

        #same vertex different times
        Q[(num_of_towns * (v1-1) + (t1+1), num_of_towns * (v2-1) + (t1+1))] += P
        Q[(num_of_towns * (v2-1) + (t1+1), num_of_towns * (v1-1) + (t1+1))] += P


# same vertex different times
#for v in G.nodes:
 #   for i, j in G.edges:
  #      Q[(num_of_towns * (i-1) + v, num_of_towns * (j-1) + v)] += P
   #     Q[(num_of_towns * (j-1) + v, num_of_towns * (i-1) + v)] += P

print(len(Q))


sampler = EmbeddingComposite(DWaveSampler(solver={'topology__type': 'chimera'}))
S = sampler.sample_qubo(Q,num_reads=numruns)


# ------- Print results to user -------

print('-' * 60)
print('{:>15s}{:>15s}{:^15s}'.format('Not Cover','Cover','Energy'))
print('-' * 60)

for sample, E in S.data(fields=['sample','energy']):
    S0 = [k for k,v in sample.items() if v == 0]        #first set (0 nodes)
    S1 = [k for k,v in sample.items() if v == 1]        #second set (1 nodes)
    print('{:>15s}{:>15s}{:^15s}'.format(str(S0),str(S1),str(E)))


# THE NEXT PART IS DISPLAYING THE BEST RESULT GRAPHICALLY
# ------- Display results to user -------
# Grab best result
# Note: "best" result is the result with the lowest energy
# Note2: the look up table (lut) is a dictionary, where the key is the node index
#   and the value is the set label. For example, lut[5] = 1, indicates that
#   node 5 is in set 1 (S1).
lut = S.first.sample

# Interpret best result in terms of nodes and edges
S0 = [node for node in G.nodes if not lut[node]]            #storing the first set
S1 = [node for node in G.nodes if lut[node]]                #storing the second set
#cut_edges = [(u, v) for u, v in G.edges if lut[u]!=lut[v]]  #storing edges that were cut (between S1 and S0)  
#uncut_edges = [(u, v) for u, v in G.edges if lut[u]==lut[v]]#storing edges that weren't cut (within each set)

# Print the solution for the user
print('order is: ', len(S1))
print(S1)

# Saving the initial graph
pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, node_size=700)
nx.draw_networkx_edges(G, pos, edgelist=G.edges, style='solid', width=3)
#nx.draw_networkx_labels(G, pos)

# node labels
nx.draw_networkx_labels(G, pos, font_size=20, font_family="sans-serif")
# edge weight labels
edge_labels = nx.get_edge_attributes(G, "weight")
nx.draw_networkx_edge_labels(G, pos, edge_labels)


filename = "Town_layout.png"
plt.savefig(filename, dpi=600)
print("\nOriginal plot has been saved to {}". format(filename))

plt.clf()   #clear old graph for new one

print(S)
