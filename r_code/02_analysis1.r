# get the subgraph of the animal-only studies (mesh) and that of human studies (mesh) 
G.color <- fix.vertex.colors(G)
G.mesh.animal <- subgraph(G.color, V(G.color)[V(G.color)$mesh.animal==TRUE])
summ.G.mesh.animal <- summarize.the.graph(G.mesh.animal, 'mesh.animal')
summ.G <- rbind(summ.G, summ.G.mesh.animal)

G.mesh.animal.only <- subgraph(G.color, V(G.color)[(V(G.color)$mesh.animal==TRUE & V(G.color)$mesh.human==FALSE)])
summ.G.mesh.animal.only <- summarize.the.graph(G.mesh.animal.only, 'mesh.animal.only')
summ.G <- rbind(summ.G, summ.G.mesh.animal.only)

G.mesh.human <- subgraph(G.color, V(G.color)[V(G.color)$mesh.human==TRUE])
summ.G.mesh.human <- summarize.the.graph(G.mesh.human, 'mesh.human')
summ.G <- rbind(summ.G, summ.G.mesh.human)


# get the subgraph of the index studies - ALL
G.all.index <- get.subgraphs.of.index.studies(G.color, which='all')
summ.G.all.index <- summarize.the.graph(G.all.index, 'all.index.studies')
summ.G <- rbind(summ.G, summ.G.all.index)

# get the subgraph of the index studies - RCTs
G.rct.index <- get.subgraphs.of.index.studies(G.color, which='rct')
summ.G.rct.index <- summarize.the.graph(G.rct.index, 'rct.index.studies')
summ.G <- rbind(summ.G, summ.G.rct.index)

# get the subgraph of the index studies - OBS
G.obs.index <- get.subgraphs.of.index.studies(G.color, which='obs')
summ.G.obs.index <- summarize.the.graph(G.obs.index, 'obs.index.studies')
summ.G <- rbind(summ.G, summ.G.obs.index)

# get the subgraph of primary data, ie. rct or obs studies based on mesh 
G.mesh.primary.data <- subgraph(G.color, V(G.color)[(V(G.color)$design=='rct' | V(G.color)$design=='obs')])
summ.G.mesh.primary.data <- summarize.the.graph(G.mesh.primary.data, 'mesh.primary.data')
summ.G <- rbind(summ.G, summ.G.mesh.primary.data)

# get the subgraph of mesh-RCTS 
G.mesh.rct <- subgraph(G.color, V(G.color)[V(G.color)$design=='rct'])
summ.G.mesh.rct <- summarize.the.graph(G.mesh.rct, 'mesh.rct')
summ.G <- rbind(summ.G, summ.G.mesh.rct)

# get the subgraph of mesh-RCTS 
G.mesh.obs <- subgraph(G.color, V(G.color)[V(G.color)$design=='obs'])
summ.G.mesh.obs <- summarize.the.graph(G.mesh.obs, 'mesh.obs')
summ.G <- rbind(summ.G, summ.G.mesh.obs)




###########
### Get the main path
###########

# this reads the pajek citation-weight vectors and attaches them 
# to the graph
spc.weights.vec <- read.pajek.vec.clu.file('spc.vec')
#splc.weights.vec <- read.pajek.vec.clu.file('splc.vec')
#spnp.weights.vec <- read.pajek.vec.clu.file('spnp.vec')

V(G)$spc.weight <- spc.weights.vec[,2]
#V(G)$splc.weight <- splc.weights.vec[,2]
#V(G)$spnp.weight <- spnp.weights.vec[,2]

# this reads the pajek main-path vector (0= not in main path, 1=in main path)
# and attached them to the graph 
spc.path.vec <- read.pajek.vec.clu.file('spc.clu')
#splc.path.vec <- read.pajek.vec.clu.file('splc.clu')
#spnp.path.vec <- read.pajek.vec.clu.file('spnp.clu')

V(G)$spc.is.path <- spc.path.vec[,2]
#V(G)$splc.is.path <- splc.path.vec[,2]
#V(G)$spnp.is.path <- spnp.path.vec[,2]


# These are the main path subgraphs directly from pajek 
# the $id attribute is the same as in G.
# Then I copy the attributes of G to the subgraphs based on the ID
# NB I can do it easy with sthg like 
# V(G)[V(G)$spc.is.path==1]$att.name = X

G.spc <- read.pajek.net.file('spc.net')
#G.splc <- read.pajek.net.file('splc.net')
#G.spnp <- read.pajek.net.file('spnp.net')

G.spc <- copy.vertex.attributes.between.graphs(G, G.spc)
#G.splc <- copy.vertex.attributes.between.graphs(G, G.splc)
#G.spnp <- copy.vertex.attributes.between.graphs(G, G.spnp)

# summarize the main path - SPC algorithm 
summ.G.spc <- summarize.the.graph(G.spc, 'SPC main path')
summ.G <- rbind(summ.G, summ.G.spc)

# write the table with the descriptives of the graphs
write.table(summ.G, file = paste(project, 'graph.summaries.txt', sep=''))


