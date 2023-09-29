# ------ Import necessary packages ----
from collections import defaultdict     #special type of python dictionary

from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import networkx as nx       #graphing package

import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt    #plotting package

# ------- Set up our graph -------

# Create empty graph
G = nx.Graph()

# Add edges to the graph (also adds nodes)
# same formulation as in the QUBO examples given in the report
G.add_edges_from([(1,2),(1,7),(2,3),(2,7),(3,4),(3,6),(4,6),(4,5),(4,7),(5,7),(6,7)])

# ------- Set up our QUBO dictionary -------

# Initialize our Q matrix
# defaultdict is a way to decalre a dictionary that let's you create labels whenever 
Q = defaultdict(int)

#un comment the next two lines if you wanna see how G.edges stores the edges (ans: tuples)
#print("The edges are: ")
#print(G.edges)

# Update Q matrix for every edge in the graph
# We try to minimize (2*x_{i}*x_{j} - x_{i} - x{j}) where i,j cycle through all the edges
# Whenever we encounter x_{i} or x{j} as standalones we add (-1) to that matrix element
# For all i-j edges we add 1 to the symmetric matrix elements
for i, j in G.edges:
    #print("{}   {}".format(i,j))
    Q[(i,i)]+= -1
    Q[(j,j)]+= -1
    Q[(i,j)]+= 1
    Q[(j,i)]+= 1

#print("The matrix after assignment is: {}".format(Q))

# ------- Run our QUBO on the QPU -------
# Set up QPU parameters
chainstrength = 8
numruns = 10

# Run the QUBO on the solver from your config file
sampler = EmbeddingComposite(DWaveSampler(solver={'topology__type': 'chimera'}))
response = sampler.sample_qubo(Q,num_reads=numruns)

# ------- Print results to user -------
print('-' * 60)
print('{:>15s}{:>15s}{:^15s}{:^15s}'.format('Set 0','Set 1','Energy','Cut Size'))
print('-' * 60)

# cycle throug the different resulting sets and their energies to display them 
for sample, E in response.data(fields=['sample','energy']):
    S0 = [k for k,v in sample.items() if v == 0]        #first set
    S1 = [k for k,v in sample.items() if v == 1]        #second set
    print('{:>15s}{:>15s}{:^15s}{:^15s}'.format(str(S0),str(S1),str(E),str(int(-1*E))))


# THE NEXT PART IS DISPLAYING THE BEST RESULT GRAPHICALLY
# ------- Display results to user -------
# Grab best result
# Note: "best" result is the result with the lowest energy
# Note2: the look up table (lut) is a dictionary, where the key is the node index
#   and the value is the set label. For example, lut[5] = 1, indicates that
#   node 5 is in set 1 (S1).
lut = response.first.sample

# Interpret best result in terms of nodes and edges
S0 = [node for node in G.nodes if not lut[node]]            #storing the first set
S1 = [node for node in G.nodes if lut[node]]                #storing the second set
cut_edges = [(u, v) for u, v in G.edges if lut[u]!=lut[v]]  #storing edges that were cut (between S1 and S0)  
uncut_edges = [(u, v) for u, v in G.edges if lut[u]==lut[v]]#storing edges that weren't cut (within each set)

# Display best result
pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, nodelist=S0, node_color='r')
nx.draw_networkx_nodes(G, pos, nodelist=S1, node_color='c')
nx.draw_networkx_edges(G, pos, edgelist=cut_edges, style='dashdot', alpha=0.5, width=3)
nx.draw_networkx_edges(G, pos, edgelist=uncut_edges, style='solid', width=3)
nx.draw_networkx_labels(G, pos)

# Saving the best result as a png
filename = "maxcut_plot.png"
plt.savefig(filename, dpi=600)
print("\nYour plot is saved to {}".format(filename))

print(response)