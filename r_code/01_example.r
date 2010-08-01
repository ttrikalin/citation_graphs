library(igraph)
source('citation_graph_functions.r')

#read the graph 
G <- read.graph('graph.xml', format='graphml')

#if it is acyclic, correct with preprints
if (is.acyclic(G)[1]==FALSE) {
    G.prime <- fix.all.preprints(G) # if it cannot be done it returns NULL 
    if (is.null(G.prime)==FALSE) {
        G <- G.prime
        G.prime<-0 # free the memory :) 
    }
}

# this prints the pajek file:
#write.pajek.mat.file(G, 'graph.pajek')

# this reads the pajek citation-weight vectors and attaches them 
# to the graph
spc <- read.pajek.vertex.weights.vec.file('spc.vec')
slpc <- read.pajek.vertex.weights.vec.file('slpc.vec')
spnp <- read.pajek.vertex.weights.vec.file('spnp.vec')

V(G)$spc.weight <- spc[,2]
V(G)$slpc.weight <- slpc[,2]
V(G)$spnp.weight <- spnp[,2]

# This gets the spc main path using a threshold t=0.01
# appending all the nodes of interest
spc.nodes <- V(G)[V(G)$spc.weight>0.01]
index.nodes <- V(G)[V(G)$is_index==1]

g <- subgraph(G, c(spc.nodes, index.nodes))
