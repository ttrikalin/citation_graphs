library(igraph)
#source('citation_graph_functions_2.r')
# take a graph and create a layout 

gimme.layout<- function(G, offset=0.1, att.name = 'color', att.val=c('red', 'blue')) {
    
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
    #add.offset <- ((get.vertex.attribute(G, att.name)==att.val) - rep(.5, length(years)))* rep(offset*Y.spread, length(years))
    # assing +.5 to those with red color
    add.offset <- (get.vertex.attribute(G, att.name)==att.val[1]) * rep(.5, length(years))
    # if contrast with blue color, assign -0.5 to blue color
    if (length(att.val)==2) {
        add.offset <- add.offset + (get.vertex.attribute(G, att.name)==att.val[2]) * rep(-.5, length(years))
    }
    add.offset <- add.offset * rep(offset*Y.spread, length(years))
    
    L <- matrix(rep(1,2*length(years)),ncol=2)
    L[,1] <- years
    L[,2] <- Y + add.offset
    
    
    return(L)
} 

is.calendar.sane <- function (G) {
    
    count <- 0
    n.edges <- ecount(G)
    M <- get.edgelist(G)
    to.delete <- NULL
    
    sane <- TRUE
    
    for (edge in 0:(n.edges-1)) {
        node1 <- M[edge+1,1]
        node2 <- M[edge+1,2]
        if (V(G)[node1]$year > V(G)[node2]$year) {
            cat("edge ", edge, " from ", node1, " to ", node2, " is not sane!\n", sep="")
            cat("year[", node1, "]=", V(G)[node1]$year, " > year [", node2, "] =", V(G)[node2]$year, "\n", sep="")
            count <- count + 1
            to.delete <- c(to.delete, edge)
            sane <- FALSE
        }
    }
    
    #print(to.delete)
    #G.corrected <- delete.edges(G, to.delete)
    
    return (c(sane=sane, violations=count, to.delete=to.delete))
}

make.calendar.sane <- function (G) {

    count <- 0
    n.edges <- ecount(G)
    M <- get.edgelist(G)
    to.delete <- NULL

    sane <- TRUE

    for (edge in 0:(n.edges-1)) {
        node1 <- M[edge+1,1]
        node2 <- M[edge+1,2]
        if (V(G)[node1]$year > V(G)[node2]$year) {
            cat("edge ", edge, " from ", node1, " to ", node2, " is not sane!\n", sep="")
            cat("year[", node1, "]=", V(G)[node1]$year, " > year [", node2, "] =", V(G)[node2]$year, "\n", sep="")
            count <- count + 1
            to.delete <- c(to.delete, edge)
            sane <- FALSE
        }
    }

    #print(to.delete)
    G.corrected <- delete.edges(G, to.delete)

    return (G.corrected)
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
    if (sum((att.list==rep('shape', length(att.list))))==0) {
        V(G)$shape <- "circle"
    }
 
   
    # This adds a preprint vertex for node1 if one does not exist 
    if (V(G)[node1]$preprint =="None") {
        G <- add.vertices(G, 1, id=paste('p', node1, sep=''), year=V(G)[node1]$year)
        V(G)[node1]$preprint=length(V(G))-1   # remember, zero based indexing
        V(G)[length(V(G))-1]$preprint="None"  # preprints have no preprints
        V(G)[length(V(G))-1]$size = V(G)[node1]$size  # get the parent size
        V(G)[length(V(G))-1]$is_index = V(G)[node1]$is_index  # get the parent index
        V(G)[length(V(G))-1]$shape = V(G)[node1]$shape  # get the parent shape
        V(G)[length(V(G))-1]$color = "NA"  # colorless
    }
    
    # This adds a preprint vertex for node2 if one does not exist 
    if (V(G)[node2]$preprint =="None") {
        G <- add.vertices(G, 1, id=paste('p', node2, sep=''), year=V(G)[node2]$year)
        V(G)[node2]$preprint=length(V(G))-1   # remember, zero based indexing
        V(G)[length(V(G))-1]$preprint="None"  # preprints have no preprints
        V(G)[length(V(G))-1]$size = V(G)[node2]$size  # get the parent size
        V(G)[length(V(G))-1]$is_index = V(G)[node2]$is_index  # get the parent index
        V(G)[length(V(G))-1]$shape = V(G)[node2]$shape  # get the parent shape
        V(G)[length(V(G))-1]$color = "NA"  # colorless
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

read.pajek.net.file <- function(filename){
    # returns a graph from the net.file
    vertex.header<-read.table(file=filename,nrows=1)
    n.vertices <- as.numeric(vertex.header[2])
    vertices<-read.table(file=filename,skip=1,nrows=n.vertices, as.is = TRUE)
    
    edges <-read.table(file=filename,skip=n.vertices+2, as.is = TRUE)
    
    # note that this graph has a node 0 which will be removed in the end
    g <- graph.edgelist(as.matrix(edges[,1:2])) 
    
    # vertex ids - one more at 0
    V(g)$id=c('dropme',vertices[,2])
    
    # edge weight - these are not more than needed as the edgelist is OK
    E(g)$weight = edges[,3]
    
    # because igraph is 0 based and the matrix started from non-zero
    g<- delete.vertices(g, 0)  
    

#    else {
#        xxrow<-read.table(file=filename,skip=1,nrows=nn[[3]],as.is=TRUE)
#        xxcol<-read.table(file=filename,skip=nn[[3]]+1, nrows=nn[[2]]-nn[[3]],fill=TRUE)
#        n<-read.table(file=filename,skip=nn[[2]]+2)
#        rownames(n)<-xxrow[[2]]
#        colnames(n)<-xxcol[[2]] 
#    }
    #as.matrix(n)
    return(g)
}


read.pajek.vec.clu.file <- function (filename) {
    
    weights <- read.table(filename, skip=1)
    vertices <- seq(0, dim(weights)[1]-1)
    W <- cbind(vertices, weights)
    colnames(W) <- c('vertex', 'weight')
    rownames(W) <- vertices
    
    return(weights=W)
    
}



parse.mesh.string <- function(mesh.string) {
    # takes a string of mesh terms separated by '#' and returns boolean 
    # on humans or animals
    
    if (is.na(mesh.string)) {
        mesh.string <- "#None#"
    }
    
    human.terms <- c('humans', 'adult', 'aged', 'middle aged','aged, 80 and over', 'male', 'female', 'adolescent', 'infant', 'infant, newborn', 'child', 'child, preschool')
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
    
    if (is.na(pub.type.string)) {
        pub.type.string <- "#None#"
    }
    
    primary.data.terms <- c('randomized controlled trial', 'clinical trial', 'cohort studies', 'prospective studies', 'longitudinal studies', 'follow-up studies', 'comparative study', 'evaluation studies')
    rct.data.terms <- c('randomized controlled trial')
    
    review.terms <- c('meta-analysis', 'metaanalysis', 'review')
    commentary.terms <- c('editorial', 'letter')
    
    primary.data <- FALSE
    rct.data <- FALSE
    review <- FALSE
    commentary <- FALSE
    
    for (pd in 1:length(primary.data.terms)) {
        if (regexpr(paste('#', primary.data.terms[pd],'#', sep=''), pub.type.string, ignore.case=TRUE)[1]!=-1) {
            primary.data <- TRUE
            break
        }
    }
    for (pd in 1:length(rct.data.terms)) {
        if (regexpr(paste('#', rct.data.terms[pd],'#', sep=''), pub.type.string, ignore.case=TRUE)[1]!=-1) {
            rct.data <- TRUE
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
    return (c(primary.data=primary.data, review=review, commentary=commentary, rct.data=rct.data))
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



copy.vertex.attributes.between.graphs <- function (fromG1, toG2) {
    # copy from G1 to G2 based on the vertex id attribute
    # G2 MUST BE A SUBGRAPH OF G1
    # not a generic function
    G1.vertices <- as.numeric(sub('n', '', V(toG2)$id))
    
    attribute.names <- list.vertex.attributes(fromG1)
    for (att in 1:length(attribute.names)) {
        toG2 <- set.vertex.attribute(toG2, attribute.names[att], V(toG2), get.vertex.attribute(fromG1, attribute.names[att], G1.vertices))        
    }
    return(toG2)
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



fix.vertex.colors <- function(G) {
    
    vertices <- V(G)
    for (i in 0:(length(vertices)-1)) {
        pub.type <- parse.pub.types.string(vertices[i]$pubtypes)
        mesh <-  parse.mesh.string(vertices[i]$mesh)
        
        if (mesh[2]==TRUE & V(G)[i]$is_index==0) { # publications on animals 
            V(G)[i]$shape='square'
            V(G)[i]$mesh.animal=TRUE
            #V(G)[i]$color='plum'
        } 
        else {
            V(G)[i]$mesh.animal=FALSE
        }
        
        if (mesh[1]==TRUE) { # publications on humans 
            V(G)[i]$mesh.human=TRUE
            #V(G)[i]$color='plum'
        }
        else {
            V(G)[i]$mesh.human=FALSE
        }
        
        if (pub.type[3]==TRUE & V(G)[i]$is_index==0) { # commentaries
            V(G)[i]$size=2
            V(G)[i]$color='black'
        }
        
        if (pub.type[2]==TRUE & V(G)[i]$is_index==0) { # reviews may overwrite commentaries
            #V(G)[i]$size=
            V(G)[i]$color='green'
        }
        if (pub.type[4]==TRUE & V(G)[i]$is_index==0 ) { # rct or obs -- will fill in the index later
            V(G)[i]$design='rct'
        }
        else if (pub.type[1]==TRUE & pub.type[4]==FALSE & V(G)[i]$is_index==0) {
            V(G)[i]$design='obs'
        }
        else if (V(G)[i]$is_index==0) {
            V(G)[i]$design='other'
        }
        
    }
    # fill in the design of the index papers
    nodes.rct <- V(G)[V(G)$color=='red']
    V(G)[nodes.rct]$design = rep('rct', length(nodes.rct))
    nodes.obs <- V(G)[V(G)$color=='blue']
    V(G)[nodes.obs]$design = rep('obs', length(nodes.obs))
    
    return(G)
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
    small.g.obs <- subgraph(G, v.list[(length(rct)+1):length(v.list)])
    
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


get.shortest.path.distribution <- function(G, nodes1=NULL, nodes2=NULL, names=c('rct', 'obs')) {
    if (is.null(nodes1)) {
        nodes1 <- V(G)[V(G)$color=='red']
    }
    if (is.null(nodes2)) {
        nodes2 <- V(G)[V(G)$color=='blue']
    }
    
    ind1 <- nodes1+1 # 0 based indexing 
    ind2 <- nodes2+1 # 0 based indexing 
    
    #print(ind1, '\n' )
    #print(ind2, '\n')
    
    M.in <- shortest.paths(G, mode='in')[ind1, ind2]
    M.out <- shortest.paths(G, mode='out')[ind1, ind2]
    M.all <- shortest.paths(G, mode='all')[ind1, ind2]
    #print(M.in)
    
    
    # return(c(M.in=M.in, M.out=M.out, M.all=M.all))
    # total paths 
    tot.paths <- length(ind1)*length(ind2)
    
    n.in.paths  <- sum(M.in!=Inf)
    n.in.bridges <- sum(M.in==1)
    n.out.paths <- sum(M.out!=Inf)
    n.out.bridges <- sum(M.out==1)
    n.all.paths <- sum(M.all!=Inf)
    n.all.bridges <- sum(M.all==1)
    
    
    if (n.in.paths >0) {
        mean.in.paths <- mean(M.in[M.in!=Inf])
        min.in.paths <- min(M.in[M.in!=Inf])
        max.in.paths <- max(M.in[M.in!=Inf])
        df.in <- c(n.in.paths, mean.in.paths, min.in.paths, max.in.paths)
        df.in <- rbind(df.in, c(n.in.bridges, n.in.bridges, n.in.bridges, n.in.bridges ))
    }
    else {
        df.in <- c(n.in.paths, NA, NA, NA)
        df.in <- rbind(df.in, c(n.in.bridges, n.in.bridges, n.in.bridges, n.in.bridges ))
    }
    
    
    if (n.out.paths >0) {
        mean.out.paths <- mean(M.out[M.out!=Inf])
        min.out.paths <- min(M.out[M.out!=Inf])
        max.out.paths <- max(M.out[M.out!=Inf])
        df.out <- c(n.out.paths, mean.out.paths, min.out.paths, max.out.paths)
        df.out <- rbind(df.out, c(n.out.bridges, n.out.bridges, n.out.bridges, n.out.bridges ))
    }
    else {
           df.out <- c(n.out.paths, NA, NA, NA)
           df.out <- rbind(df.out, c(n.out.bridges, n.out.bridges, n.out.bridges, n.out.bridges ))
    }
    
    
    if (n.all.paths >0 ) {
        mean.all.paths <- mean(M.all[M.all!=Inf])
        min.all.paths <- min(M.all[M.all!=Inf])
        max.all.paths <- max(M.all[M.all!=Inf])
        df.all <- c(n.all.paths, mean.all.paths, min.all.paths, max.all.paths)
        df.all <- rbind(df.all, c(n.all.bridges, n.all.bridges, n.all.bridges, n.all.bridges ))
    }
    else {
        df.all <- c(n.all.paths, NA, NA, NA)
        df.all <- rbind(df.all, c(n.all.bridges, n.all.bridges, n.all.bridges, n.all.bridges ))
    }
    
    df <- rbind(c(tot.paths, NA, NA, NA), df.in, df.out, df.all)
    colnames(df) <- c('number', 'mean', 'min', 'max')
    rownames(df) <- c('possible.shortest.paths', paste(names[1],'.in.paths', sep=''), paste(names[1],'.in.bridges', sep=''), paste(names[2],'.in.paths', sep=''), paste(names[2],'.in.bridges', sep=''), 'undirected.paths', 'undirected.bridges')
    return(df)
}

summarize.the.graph <- function (G, topicname='a graph') {
    
    n.v <- vcount(G)
    n.e <- ecount(G)
    d.G <- graph.density(G)
    #n.rct <- length(clean.my.python.list.string(get.graph.attribute(G,'rct_pmid_list')))
    #n.obs <- length(clean.my.python.list.string(get.graph.attribute(G,'obs_pmid_list')))
    n.rct <- sum(V(G)$color=='red')
    n.obs <- sum(V(G)$color=='blue')
    
    G.color <- fix.vertex.colors(G)
    n.mesh.animal <- sum(V(G.color)$shape=='square')
    n.mesh.review <- sum(V(G.color)$color=='green')
    years <- summary(V(G)$year)[c(1,3,5)]
    
    m <- data.frame(t(c(topicname, n.v, n.e, d.G, n.rct, n.obs, n.mesh.animal, n.mesh.review, years )), stringsAsFactors=FALSE)
    colnames(m) <- c('topic name','n.vertices', 'n.edges', 'density', 'index.rcts', 'index.obs', 'n.mesh.animal', 'n.mesh.review', 'years.min', 'years.median', 'years.max')
    return(m)
}
 

myplot <- function (G, step=3, offset=0, curved=.1, width=.5, v.size=NULL) {
    
    L <- gimme.layout(G, offset=offset)
    lmin <- min(L[,1])
    lmax <- max(L[,1])
    if (is.null(v.size)){
        v.size <- V(G)$size
    }
    
    at <- seq(from=-1, to=1, by=2*step/(lmax-lmin))
    labels <- seq(from=lmin, length.out=(lmax-lmin)/step, by=step)
    if (length(labels)<length(at)) {
        labels <- seq(from=lmin, length.out=(lmax-lmin)/step+1, by=step)
    }
    plot(G, layout=L, edge.curved=curved, edge.arrow.size=0.2, edge.width=width, vertex.size=v.size)
    axis(1, at=at, labels=labels)
}

