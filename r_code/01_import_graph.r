# start with a clean slate 
#LS1 <- ls()
#LS2 <- NULL
#for (i in 1:length(LS1)) {
#    if (LS1[i] != 'cdf') {
#        LS2 <- c(LS2, LS1[i])
#    }
#    LS2 <- c(LS2, 'LS1')
#}
#remove(list=LS2)
#remove(LS2)


# load functions 
library(igraph)
source('citation_graph_functions.r')

# set project title
project <- 'vitE_cvd'


# form the GMLpath
XMLpath <- paste('../Results/', project,'/graph', sep='')

#read the graph 
oldwd <- getwd()
setwd(XMLpath)
G <- read.graph('graph.xml', format='graphml')

# make sure circle is the default shape 
V(G)$shape <- 'circle'

#if it is acyclic, correct with preprints
if (is.acyclic(G)[1]==FALSE) {
    G.prime <- fix.all.preprints(G) # if it cannot be done it returns NULL 
    if (is.null(G.prime)==FALSE) {
        G <- G.prime
        G.prime<-0 # free the memory :) 
    }
}

#make calendar sane
G <- make.calendar.sane(G)

# summarize the graph
summ.G <- summarize.the.graph(G, 'main')

# this prints the pajek file:
write.pajek.mat.file(G, 'graph.pajek')

# go the the analysis part  - descriptives of graphs 
source('../../../r_code/02_analysis1.r')

# this is some graphical output 
source('../../../r_code/03_analysis.r')

# this reports some communication statistics 
source('../../../r_code/04_shortest_paths.r')

