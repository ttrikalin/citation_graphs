# this gets the subgraph made of index vertices
#

extract.subgraph <- function(G, index.vertices, mode='all') {
	g.vertices <- NULL
	n <- length(index.vertices) 
	for (v in 0:(n-1)) {
		g.vertices <- c(g.vertices, subcomponent(G, v, mode))
	}

	g <- subgraph(G, g.vertices)
	return(g)
}


