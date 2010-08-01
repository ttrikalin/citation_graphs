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




is.acyclic <- function(G) {
    order <- diameter(G)
    acyclic <- TRUE
    for (v in 0:(length(V(G))-1)) {
        v.in <- neighborhood(G, order, v, mode='in')[[1]]
        v.out <- neighborhood(G, order, v, mode='out')[[1]]
        n.in <- length(v.in)
        n.out <- length(v.out)
        
        if ((n.in>1) & (n.out>1)) {       # if only the index node in the 'in' or 'out', no problem
            for (v1 in 2:length(v.in)) {  # the first element is always the index node
                for (v2 in 2:length(v.out)) {
                    if ((v.in[v1]==v.out[v2])) {
                        cat("Found cyclic path for node", v, '!\n')
                        cat("incoming:", v.in[2:n.in], '\n')
                        cat("outgoing:", v.out[2:n.out], '\n')
                        acyclic <- FALSE
                        return(c(acyclic=acyclic, node1=v, node2=v.in[v1]))
                    }
                }
            }            
        }   
    }
    return(acyclic=acyclic)
}

fix.all.preprints <- function (G) {
    goon <- TRUE
    count <- 0 

    G<- simplify(G)   # just in case

    while (goon == TRUE) {
    
        print(count)
        
        check <- is.acyclic(G)
    
        if (check[1]==TRUE) {
            if (count == 0) {
                print("Graph was acyclic, no modification")
            }
            return(G)
        }
        
        G.prime <- fix.mutual(G, check[2], check[3])
        if (is.null(G.prime)) {
            cat("Error when checking nodes ", check[2], " and ", check[3], '\n', sep='' )
            return(NULL)
        }
        G <- G.prime
        count <- count + 1
    }
    return(G)
}

fix.mutual <- function(G, node1, node2) {
    # fixes the preprint problem for a set of mutual nodes
    
    # if the nodes are not mutual, abort
    if (is.mutual(G, E(G, c(node1, node2)))==FALSE) {
        cat('Nodes are not mutual.\nfix.mutual(G, ', node1, ', ', node2, ') cannot fix the cyclic path', sep='')
        return(NULL)
    }
    
    # just in case, simplify
    G <- simplify(G)   
    
    #  first check if there is a preprint or an attribute
    att.list <- list.vertex.attributes(G)
    if (sum((att.list==rep('preprint', length(att.list))))==0) {
        V(G)$preprint <- "None"
    }
 
   
    # This adds a preprint vertex for node1 if one does not exist 
    if (V(G)[node1]$preprint =="None") {
        G <- add.vertices(G, 1, id=paste('p', node1, sep=''), year=V(G)[node1]$year)
        V(G)[node1]$preprint=length(V(G))-1   # remember, zero based indexing
        V(G)[length(V(G))-1]$preprint="None"  # preprints have no preprints
        V(G)[length(V(G))-1]$size = V(G)[node1]$size  # get the parent size
        V(G)[length(V(G))-1]$is_index = V(G)[node1]$is_index  # get the parent index
    }
    
    # This adds a preprint vertex for node2 if one does not exist 
    if (V(G)[node2]$preprint =="None") {
        G <- add.vertices(G, 1, id=paste('p', node2, sep=''), year=V(G)[node2]$year)
        V(G)[node2]$preprint=length(V(G))-1   # remember, zero based indexing
        V(G)[length(V(G))-1]$preprint="None"  # preprints have no preprints
        V(G)[length(V(G))-1]$size = V(G)[node2]$size  # get the parent size
        V(G)[length(V(G))-1]$is_index = V(G)[node2]$is_index  # get the parent index
    }
     
    # This makes the following
    # Transform:
    #  A -> B
    #  B -> A
    #
    # To:
    #  A -> pB
    #  B -> pA
    #  A -> pA
    #  B -> pB
    
    G <- delete.edges(G, E(G, c(node1, node2)))
    G <- delete.edges(G, E(G, c(node2, node1)))
    G <- add.edges(G, c(node1, as.numeric(V(G)[node1]$preprint)))
    G <- add.edges(G, c(node1, as.numeric(V(G)[node2]$preprint)))
    G <- add.edges(G, c(node2, as.numeric(V(G)[node1]$preprint)))
    G <- add.edges(G, c(node2, as.numeric(V(G)[node2]$preprint)))
    
    # if any of the nodes had a preprint from a previous correction you'll have 2 edges
    # simplify!
    G <- simplify(G)
    
    return(G)
}


write.pajek.mat.file <- function(G, filename , name='id', VERBOSE=FALSE) {
    
    rs <- paste('*Vertices ', length(V(G)), '\r\n', sep='')
    for (v in 0:(length(V(G))-1)) {
        rs <- paste(rs, v+1, ' "', get.vertex.attribute(G, name, v), '"\r\n', sep='')
    }
    rs <- paste(rs, '*Arcs ', length(E(G)),'\r\n', sep='') 
    
    for (e in 0:(length(E(G))-1)) {   
        rs <- paste(rs, get.edges(G, E(G)[e])[1]+1, ' ',get.edges(G, E(G)[e])[2]+1, '\r\n',sep='')
    }
    
    write(rs, file=paste(filename, '.mat', sep=''))
    
    if (VERBOSE==TRUE) {
        return(rs)
    }
    else {
        return("Success!")
    }
}

read.pajek.vertex.weights.vec.file <- function (filename) {
    
    weights <- read.table(filename, skip=1)
    vertices <- seq(0, length(vector)-1)
    W <- cbind(vertices, weights)
    colnames(W) <- c('vertex', 'weight')
    
    return(weights=W)
    
}

keep.edges <- function (G, edge.sequence){
    
    G <- delete.edges(G, E(G))
    print("here")
    G <- add.edge(G, edge.sequence)
    return(G)
}

parse.mesh.string <- function(mesh.string) {
    # takes a string of mesh terms separated by '#' and returns boolean 
    # on humans or animals
    
    human.terms <- c('humans', 'adult', 'aged', 'middle aged', 'male', 'female', 'adolescent', 'infant', 'infant, newborn', 'child', 'child, preschool')
    animal.terms <- c('animals')
    
    humans <- FALSE
    animals <- FALSE
    for (ht in 1:length(human.terms)) {
        if (regexpr(paste('#', human.terms[ht],'#', sep=''), mesh.string, ignore.case=TRUE)[1]!=-1) {
            humans <- TRUE
            break
        }
    }
    for (at in 1:length(animal.terms)) {
        if (regexpr(paste('#', animal.terms[at],'#', sep=''), mesh.string, ignore.case=TRUE)[1]!=-1) {
            animals <- TRUE
            break
        }
    }
    return (c(humans=humans, animals=animals))
}


parse.pub.types.string <- function(pub.type.string) {
    # takes a string of mesh terms separated by '#' and returns boolean 
    # on humans or animals
    
    primary.data.terms <- c('randomized controlled trial', 'clinical trial')
    review.terms <- c('meta-analysis', 'metaanalysis', 'review')
    commentary.terms <- c('editorial', 'letter')
    
    primary.data <- FALSE
    review <- FALSE
    commentary <- FALSE
    
    for (pd in 1:length(primary.data.terms)) {
        if (regexpr(paste('#', primary.data.terms[pd],'#', sep=''), pub.type.string, ignore.case=TRUE)[1]!=-1) {
            primary.data <- TRUE
            break
        }
    }
    for (rt in 1:length(review.terms)) {
        if (regexpr(paste('#', review.terms[rt],'#', sep=''), pub.type.string, ignore.case=TRUE)[1]!=-1) {
            review <- TRUE
            break
        }
    }
    for (ct in 1:length(commentary.terms)) {
        if (regexpr(paste('#', commentary.terms[ct],'#', sep=''), pub.type.string, ignore.case=TRUE)[1]!=-1) {
            commentary <- TRUE
            break
        }
    }
    return (c(primary.data=primary.data, review=review, commentary=commentary))
}

write.pajek.file.windows <- function(G,filename,name='pmid',twomode=1){
    
    # this will work only on windows machines, because the EOL in other machines
    # is \n and not \r\n as pajek expects
    
    M <- get.adjacency(G)
    if (is.null(name)==FALSE) {
        rownames(M) <- get.vertex.attribute(G, name)
    }
    else {
        rownames(M) <- paste('v',seq(0, dim(M)[1]-1), sep='')
    }
    if ((dim(M)[1] == dim(M)[2]) & (twomode!=2)) {
        write(paste("*Vertices",dim(M)[1]), file = filename);
        write(paste(seq(1,length=dim(M)[1]),' "',rownames(M),
            '"',sep=""), file = filename,append=TRUE);
        write("*Arcs", file = filename,append=TRUE);
        for (i in 1:dim(M)[1]) {
            for (j in 1:dim(M)[2]) {
                if (M[i,j]!=0) {
                    write(paste(i,j,M[i,j]),
                        file = filename,append=TRUE)
                }
            }
        }
    } 
#    else {
#        write(paste("*Vertices",sum(dim(M)),dim(M)[1]),
#            file = filename);
#        write(paste(1:dim(M)[1],' "',rownames(M),'"',sep=""),
#            file = filename,append=TRUE);
#        write(paste(seq(dim(M)[1]+1,length=dim(M)[2]),' "',
#            colnames(M),'"',sep=""), file = filename,append=TRUE);
#        write("*Edges", file = filename, append=TRUE);
#        for (i in 1:dim(M)[1]) {
#            for (j in 1:dim(M)[2]) {
#                if (M[i,j]!=0) {
#                    write(paste(i,j+dim(M)[1],M[i,j]),
#                        file = filename,append=TRUE)
#                }
#            }
#        }
#    }
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

