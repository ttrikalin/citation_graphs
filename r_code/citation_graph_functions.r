library(igraph)

# take a graph and create a layout 

gimme.layout<- function(G, offset=0.1, att.name = 'color', att.val='red') {
    
    years <- get.vertex.attribute(G,'year')
    
    min.year <- min(years)
    max.year <- max(years)
    how.many.per.year <- NULL
    for (i in min.year:max.year) {
        how.many.per.year <- c(how.many.per.year,sum(years==i))
    }

    Y.spread<-max(how.many.per.year)
    Y <- NULL
    for (i in 1:length(how.many.per.year)) {
        if (how.many.per.year[i]>0) {
            y <- seq.int(1,how.many.per.year[i], 1)  
            Y <- c(Y, (y-mean(y) + runif(1)/5)/Y.spread)
        }
    }
    
    
    for (v in 0:(length(years)-1)) {
        vertex.year<-get.vertex.attribute(G,'year',v)
        number.this.year <- how.many.per.year[vertex.year - min.year + 1]
        
    }
    
    #moves rcts up and obs down by 0.5*offset*Y.spread
    add.offset <- ((get.vertex.attribute(g, att.name)==att.val) - rep(.5, length(years)))* rep(offset*Y.spread, length(years))
    
    L <- matrix(rep(1,2*length(years)),ncol=2)
    L[,1] <- years
    L[,2] <- Y + add.offset
    
    
    return(L)
} 


get.vertex.index.from.attribute <- function(G, attribute.value, attribute.name = 'pmid'){
    vertices <- V(G)
    for (i in 0:(length(vertices)-1)){
        att <- get.vertex.attribute(G, attribute.name, i)
        if ((att == attribute.value)){
            return(i)
        }
    }
    
    cat("Found no vertex with the specified attribute\n")
    return(NULL)
}

clean.my.python.list.string <- function (S) {
    S <- sub("'\\]","", S)
    S <- sub("\\['","", S)
    S.list <- strsplit(S,"', '")
    return(S.list[[1]])
}

get.subgraphs.of.index.studies <- function(G, which='all') {
    
    rct <- clean.my.python.list.string(G$rct_pmid_list)
    obs <- clean.my.python.list.string(G$obs_pmid_list)    
    all.studies <- c(rct, obs)

    v.list <- NULL
    for (i in 1:length(all.studies)){
        v.list[i] <- get.vertex.index.from.attribute(G,all.studies[i])
    }
    
    small.g <- subgraph(G, v.list)
    small.g.rct <- subgraph(G, v.list[1:length(rct)])
    small.g.obs <- subgraph(G, v.list[length(rct):length(v.list)])
    
    if (which == 'all'){
        rs <- small.g
    }
    else if (which == 'rct'){
        rs <- small.g.rct
    } 
    else if (which == 'obs'){
        rs <- small.g.obs
    } 

    return(rs)
}

