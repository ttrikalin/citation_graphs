library(igraph)
source('citation_graph_functions.r')

#read the graph 
G <- read.graph('graph.xml', format='graphml')

#read the years 
years <- get.vertex.attribute(G, 'year')

#create layout 
#L <- layout.spring(G)

# years in the x-axis
#L2 <- L
#L2[,1]<-years  

# fix the plot
plot(G, edge.arrow.size=0.5, xlab='year of publication')

 